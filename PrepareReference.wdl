version 1.0

## PrepareReference.wdl
##
## Builds the GLIMPSE2 reference panel for low-pass imputation, one shard per
## chromosome. Reproduces the chr22 pilot procedure that yielded 167,215 sites.
##
## Site selection is the union of two criteria on the HGDP+1KGP panel:
##   - MAF >= 0.01 across all 4,151 samples
##   - MAF >= 0.01 within the 808 East Asian samples
## The union matters because roughly 13% of EAS-common variants sit below 1%
## globally and would otherwise vanish from a Korean cohort.
##
## Per chromosome this emits sites.vcf.gz (drives chunking and mpileup
## targeting), sites.tsv.gz (for bcftools call -T), panel.bcf (genotypes at
## those sites), chunks.txt, and a tarball of the binary reference files that
## GLIMPSE2_phase consumes.

workflow PrepareReference {
  input {
    Array[String] chroms = ["chr22"]

    String panel_prefix = "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp"
    File   eas_samples
    String gmap_dir

    Float  maf_min = 0.01

    String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    String glimpse_docker  = "simrub/glimpse:v2.0.0-27-g0919952_20221207"
  }

  scatter (c in chroms) {
    call ExtractSites {
      input:
        chrom           = c,
        panel_vcf       = panel_prefix + "." + c + ".vcf.bgz",
        panel_vcf_index = panel_prefix + "." + c + ".vcf.bgz.tbi",
        eas_samples     = eas_samples,
        maf_min         = maf_min,
        docker          = bcftools_docker
    }

    call ChunkAndSplit {
      input:
        chrom       = c,
        sites_vcf   = ExtractSites.sites_vcf,
        sites_index = ExtractSites.sites_vcf_index,
        panel_bcf   = ExtractSites.panel_bcf,
        panel_index = ExtractSites.panel_bcf_index,
        gmap        = gmap_dir + "/" + c + ".b38.gmap.gz",
        docker      = glimpse_docker
    }
  }

  output {
    Array[File] out_sites_vcf       = ExtractSites.sites_vcf
    Array[File] out_sites_vcf_index = ExtractSites.sites_vcf_index
    Array[File] out_sites_tsv       = ExtractSites.sites_tsv
    Array[File] out_sites_tsv_index = ExtractSites.sites_tsv_index
    Array[File] out_panel_bcf       = ExtractSites.panel_bcf
    Array[File] out_panel_bcf_index = ExtractSites.panel_bcf_index
    Array[Int]  out_n_sites         = ExtractSites.n_sites
    Array[File] out_chunks          = ChunkAndSplit.chunks
    Array[File] out_ref_bins        = ChunkAndSplit.ref_bins
  }
}


task ExtractSites {
  input {
    String chrom
    File   panel_vcf
    File   panel_vcf_index
    File   eas_samples
    Float  maf_min
    String docker
  }

  # Arithmetic stays out of the command placeholders: the WDL parser reads
  # "1.0" inside ~{ } as member access on the integer 1 and rejects the file.
  Float maf_max = 1.0 - maf_min
  Int   disk_gb = ceil(size(panel_vcf, "GB") * 2 + 60)

  command <<<
    set -euo pipefail

    # bcftools +fill-tags derives a per-group AF when given a sample-to-group
    # table. Samples absent from the table contribute nothing to AF_eas.
    awk '{print $1"\teas"}' ~{eas_samples} > eas_group.txt
    wc -l eas_group.txt

    # Single pass over the panel: biallelic SNPs, add AF_eas, then keep sites
    # common in either the whole panel or the East Asian subset.
    bcftools view -m2 -M2 -v snps -Ou ~{panel_vcf} \
      | bcftools +fill-tags -Ou -- -t AF -S eas_group.txt \
      | bcftools view \
          -i 'INFO/AF >= ~{maf_min} && INFO/AF <= ~{maf_max} || INFO/AF_eas >= ~{maf_min} && INFO/AF_eas <= ~{maf_max}' \
          -Ob -o filtered.bcf
    bcftools index -f filtered.bcf

    N=$(bcftools index -n filtered.bcf)
    echo "$N" > n_sites.txt
    echo "sites kept on ~{chrom}: $N"

    # Genotypes only. The gnomAD annotations would otherwise drag ~55 GB of
    # dead weight through every downstream phase shard.
    bcftools annotate -x INFO,^FORMAT/GT -Ob -o panel_~{chrom}.bcf filtered.bcf
    bcftools index -f panel_~{chrom}.bcf

    bcftools view -G -Oz -o sites_~{chrom}.vcf.gz filtered.bcf
    bcftools index -f -t sites_~{chrom}.vcf.gz

    # bcftools call -C alleles -T wants REF and comma-joined ALT in one column.
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' filtered.bcf \
      | bgzip -c > sites_~{chrom}.tsv.gz
    tabix -s1 -b2 -e2 -f sites_~{chrom}.tsv.gz

    ls -lh panel_~{chrom}.bcf sites_~{chrom}.vcf.gz sites_~{chrom}.tsv.gz
  >>>

  runtime {
    docker: docker
    cpu: 8
    memory: "32 GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 0
    maxRetries: 1
  }

  output {
    File sites_vcf       = "sites_~{chrom}.vcf.gz"
    File sites_vcf_index = "sites_~{chrom}.vcf.gz.tbi"
    File sites_tsv       = "sites_~{chrom}.tsv.gz"
    File sites_tsv_index = "sites_~{chrom}.tsv.gz.tbi"
    File panel_bcf       = "panel_~{chrom}.bcf"
    File panel_bcf_index = "panel_~{chrom}.bcf.csi"
    Int  n_sites         = read_int("n_sites.txt")
  }
}


task ChunkAndSplit {
  input {
    String chrom
    File   sites_vcf
    File   sites_index
    File   panel_bcf
    File   panel_index
    File   gmap
    String docker
  }

  Int disk_gb = ceil(size(panel_bcf, "GB") * 4 + 50)

  command <<<
    set -euo pipefail

    # htslib looks for an index beside its data file, so link both together.
    ln -s ~{sites_vcf}   ./sites.vcf.gz
    ln -s ~{sites_index} ./sites.vcf.gz.tbi
    ln -s ~{panel_bcf}   ./panel.bcf
    ln -s ~{panel_index} ./panel.bcf.csi

    GLIMPSE2_chunk_static \
      --input ./sites.vcf.gz \
      --region ~{chrom} \
      --map ~{gmap} \
      --sequential \
      --threads 4 \
      --output chunks_~{chrom}.txt

    echo "chunks:"
    cat chunks_~{chrom}.txt

    # Column 3 is the buffered region GLIMPSE2 reads, column 4 the region it
    # will emit. split_reference needs both.
    while IFS=$'\t' read -r ID CHR IRG ORG REST; do
      GLIMPSE2_split_reference_static \
        --reference ./panel.bcf \
        --map ~{gmap} \
        --input-region  "$IRG" \
        --output-region "$ORG" \
        --threads 4 \
        --output reference_~{chrom}
    done < chunks_~{chrom}.txt

    ls -lh reference_~{chrom}*.bin
    ls reference_~{chrom}*.bin | wc -l

    # One tarball instead of a glob: nested File arrays trip some parsers, and
    # phase will untar these anyway.
    mkdir -p refbins
    mv reference_~{chrom}*.bin refbins/
    tar czf refbins_~{chrom}.tar.gz refbins
  >>>

  runtime {
    docker: docker
    cpu: 8
    memory: "32 GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 0
    maxRetries: 1
  }

  output {
    File chunks   = "chunks_~{chrom}.txt"
    File ref_bins = "refbins_~{chrom}.tar.gz"
  }
}
