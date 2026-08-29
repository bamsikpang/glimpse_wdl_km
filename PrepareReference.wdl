version 1.0

## PrepareReference.wdl
##
## Builds the GLIMPSE2 reference panel for low-pass imputation, one shard per
## chromosome, from the gnomAD v3.1.2 HGDP+1KGP release (4,151 samples).
##
## Site selection keeps biallelic SNPs that pass gnomAD's own filters and are
## common in either stratum:
##   - MAF >= 0.01 across all 4,151 samples, or
##   - MAF >= 0.01 within the 808 East Asian samples
## The union matters because roughly 13% of EAS-common variants sit below 1%
## globally and would otherwise vanish from a Korean cohort.
##
## Both GLIMPSE2_chunk and GLIMPSE2_split_reference abort outright unless
## AC and AN are present in INFO, so both the sites VCF and the panel BCF
## carry exactly those two fields and nothing else. Stripping gnomAD's own
## annotations before recomputing keeps the sites VCF at a few megabytes
## instead of well over a gigabyte.

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

    # A bare "command not found" surfaces as exit 127 with no hint as to
    # which tool was absent, so name them up front.
    for t in bcftools bgzip tabix; do
      command -v $t || echo "MISSING: $t"
    done
    bcftools plugin -l 2>/dev/null | grep -i fill-tags || echo "MISSING: fill-tags plugin"

    # bcftools +fill-tags derives a per-group AF from a sample-to-group table.
    # Samples absent from the table contribute nothing to AF_eas.
    awk '{print $1"\teas"}' ~{eas_samples} > eas_group.txt
    wc -l eas_group.txt

    # One pass over the panel, and the order of these three steps matters.
    # -f PASS drops what gnomAD flagged (AC0, AS_VQSR, InbreedingCoeff),
    # which otherwise pads the site count with low-quality calls near the
    # centromere. annotate -x INFO clears gnomAD's per-population fields so
    # that fill-tags recomputes AF and AF_eas from scratch and the filter
    # below reads freshly derived values.
    bcftools view -m2 -M2 -v snps -f PASS -Ou ~{panel_vcf} \
      | bcftools annotate -x INFO -Ou \
      | bcftools +fill-tags -Ou -- -t AF -S eas_group.txt \
      | bcftools view \
          -i 'INFO/AF >= ~{maf_min} && INFO/AF <= ~{maf_max} || INFO/AF_eas >= ~{maf_min} && INFO/AF_eas <= ~{maf_max}' \
          -Ob -o filtered.bcf
    bcftools index -f filtered.bcf

    N=$(bcftools index -n filtered.bcf)
    echo "$N" > n_sites.txt
    echo "sites kept on ~{chrom}: $N"

    # Strip INFO first, then let fill-tags derive AC/AN from exactly the
    # genotypes that remain. Computing them before the strip left the two
    # out of step, and split_reference cross-checks AC/AN against GT.
    # Missing calls are set to reference: GLIMPSE2 expects a complete
    # haplotype panel, and any './.' left in place makes AN disagree.
    bcftools annotate -x INFO,^FORMAT/GT -Ou filtered.bcf \
      | bcftools +setGT -Ou -- -t . -n 0p \
      | bcftools +fill-tags -Ob -o panel_~{chrom}.bcf -- -t AC,AN
    bcftools index -f panel_~{chrom}.bcf

    # The sites VCF is the panel with genotypes dropped, so its AC/AN match
    # the panel's by construction.
    bcftools view -G -Oz -o sites_~{chrom}.vcf.gz panel_~{chrom}.bcf
    bcftools index -f -t sites_~{chrom}.vcf.gz

    # bcftools call -C alleles -T wants REF and comma-joined ALT in one column.
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' filtered.bcf \
      | bgzip -c > sites_~{chrom}.tsv.gz
    tabix -s1 -b2 -e2 -f sites_~{chrom}.tsv.gz

    echo "=== INFO check ==="
    bcftools view -h sites_~{chrom}.vcf.gz | grep "^##INFO" || true
    bcftools view -h panel_~{chrom}.bcf    | grep "^##INFO" || true

    ls -lh panel_~{chrom}.bcf sites_~{chrom}.vcf.gz sites_~{chrom}.tsv.gz
  >>>

  runtime {
    docker: docker
    cpu: 8
    memory: "32 GB"
    disks: "local-disk ~{disk_gb} SSD"
    # Streaming 55 GB through piped bcftools stages runs ~1.5 h, long enough
    # that preemptible VMs were reclaimed twice before finishing.
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

    # Release tarballs name the binaries GLIMPSE2_chunk_static; images built
    # from source drop the suffix. Resolve at runtime and print what turned
    # up so a miss is diagnosable from the log.
    ls -1 /usr/local/bin /usr/bin /opt/*/bin 2>/dev/null | grep -i glimpse || true
    CHUNK=$(command -v GLIMPSE2_chunk_static || command -v GLIMPSE2_chunk)
    SPLIT=$(command -v GLIMPSE2_split_reference_static || command -v GLIMPSE2_split_reference)
    echo "chunk=$CHUNK  split=$SPLIT"

    # htslib looks for an index beside its data file, so link both together.
    ln -s ~{sites_vcf}   ./sites.vcf.gz
    ln -s ~{sites_index} ./sites.vcf.gz.tbi
    ln -s ~{panel_bcf}   ./panel.bcf
    ln -s ~{panel_index} ./panel.bcf.csi

    "$CHUNK" \
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
      "$SPLIT" \
        --reference ./panel.bcf \
        --map ~{gmap} \
        --input-region  "$IRG" \
        --output-region "$ORG" \
        --threads 4 \
        --output reference_~{chrom}
    done < chunks_~{chrom}.txt

    ls -lh reference_~{chrom}*.bin
    ls reference_~{chrom}*.bin | wc -l

    # One tarball instead of a glob: nested File arrays trip some parsers,
    # and phase unpacks these anyway.
    mkdir -p refbins
    mv reference_~{chrom}*.bin refbins/
    tar czf refbins_~{chrom}.tar.gz refbins
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "32 GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    File chunks   = "chunks_~{chrom}.txt"
    File ref_bins = "refbins_~{chrom}.tar.gz"
  }
}
