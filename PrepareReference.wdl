version 1.0

## PrepareReference.wdl
##
## Builds the GLIMPSE2 reference panel for low-pass imputation, one shard per
## chromosome. Reproduces the chr22 pilot procedure that yielded 167,215 sites.
##
## Site selection is the union of two criteria on the HGDP+1KGP panel:
##   - MAF >= 0.01 in the full 4,151-sample panel
##   - MAF >= 0.01 in the 808 East Asian samples
## The union matters because ~13% of EAS-common variants fall below 1% globally
## and would otherwise be dropped from a Korean cohort.
##
## Outputs per chromosome: sites.vcf.gz (for chunking and mpileup targeting),
## sites.tsv.gz (for bcftools call -T), panel.bcf (genotypes at those sites),
## chunks.txt, and the binary reference files consumed by GLIMPSE2_phase.

workflow PrepareReference {
  input {
    Array[String] chroms = ["chr22"]

    String panel_prefix = "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp"
    File   eas_samples
    String gmap_dir

    Float  maf_min = 0.01

    String bcftools_docker = "staphb/bcftools:1.19"
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
    Array[File]       sites_vcf   = ExtractSites.sites_vcf
    Array[File]       sites_index = ExtractSites.sites_vcf_index
    Array[File]       sites_tsv   = ExtractSites.sites_tsv
    Array[File]       sites_tsv_index = ExtractSites.sites_tsv_index
    Array[File]       panel_bcf   = ExtractSites.panel_bcf
    Array[File]       panel_index = ExtractSites.panel_bcf_index
    Array[Int]        n_sites     = ExtractSites.n_sites
    Array[File]       chunks      = ChunkAndSplit.chunks
    Array[Array[File]] ref_bins   = ChunkAndSplit.ref_bins
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

  Float  vcf_gb  = size(panel_vcf, "GB")
  Int    disk_gb = ceil(vcf_gb * 2 + 60)

  command <<<
    set -euo pipefail

    # bcftools +fill-tags computes per-group AF when handed a two-column
    # sample-to-group file. Everyone not listed is ignored for the group tag,
    # so only the EAS samples contribute to AF_eas.
    awk '{print $1"\teas"}' ~{eas_samples} > eas_group.txt
    wc -l eas_group.txt

    # One pass over the panel: restrict to biallelic SNPs, add AF_eas, then
    # keep any site common in either the full panel or the EAS subset.
    bcftools view -m2 -M2 -v snps -Ou ~{panel_vcf} \
      | bcftools +fill-tags -Ou -- -t AF -S eas_group.txt \
      | bcftools view \
          -i 'INFO/AF >= ~{maf_min} && INFO/AF <= ~{1.0 - maf_min} || INFO/AF_eas >= ~{maf_min} && INFO/AF_eas <= ~{1.0 - maf_min}' \
          -Ob -o filtered.bcf
    bcftools index -f filtered.bcf

    N=$(bcftools index -n filtered.bcf)
    echo "$N" > n_sites.txt
    echo "sites kept on ~{chrom}: $N"

    # Genotypes only. Everything downstream reads GT; the gnomAD annotations
    # would otherwise carry ~55 GB of dead weight through every phase shard.
    bcftools annotate -x INFO,^FORMAT/GT -Ob -o panel_~{chrom}.bcf filtered.bcf
    bcftools index -f panel_~{chrom}.bcf

    # Sites-only VCF drives GLIMPSE2_chunk and targets bcftools mpileup.
    bcftools view -G -Oz -o sites_~{chrom}.vcf.gz filtered.bcf
    bcftools index -f -t sites_~{chrom}.vcf.gz

    # The TSV form is what `bcftools call -C alleles -T` expects: REF and the
    # comma-joined ALT in a single column.
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' filtered.bcf \
      | bgzip -c > sites_~{chrom}.tsv.gz
    tabix -s1 -b2 -e2 -f sites_~{chrom}.tsv.gz

    ls -lh panel_~{chrom}.bcf sites_~{chrom}.vcf.gz sites_~{chrom}.tsv.gz
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    File sites_vcf       = "sites_" + chrom + ".vcf.gz"
    File sites_vcf_index = "sites_" + chrom + ".vcf.gz.tbi"
    File sites_tsv       = "sites_" + chrom + ".tsv.gz"
    File sites_tsv_index = "sites_" + chrom + ".tsv.gz.tbi"
    File panel_bcf       = "panel_" + chrom + ".bcf"
    File panel_bcf_index = "panel_" + chrom + ".bcf.csi"
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

    # Index files must sit beside their data files for htslib to find them.
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

    echo "chunks:"; cat chunks_~{chrom}.txt

    # Column 3 is the buffered region, column 4 the region GLIMPSE2 will
    # actually emit. split_reference wants both.
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
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "32 GB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    File        chunks   = "chunks_" + chrom + ".txt"
    Array[File] ref_bins = glob("reference_" + chrom + "*.bin")
  }
}
