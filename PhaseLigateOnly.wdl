version 1.0

## PhaseLigateOnly.wdl
##
## Picks up from an already-merged GL BCF and runs the two remaining steps:
## GLIMPSE2_phase over each chunk, then GLIMPSE2_ligate to stitch them back
## into a whole-chromosome callset.
##
## Splitting this out of PhaseLigate lets a merge that already succeeded be
## reused rather than repeating the ~2h40m of localizing 4,840 per-sample
## files.
##
## Two things about the environment shape the shell below. The bcftools image
## ships BusyBox coreutils, and `set -o pipefail` turns a SIGPIPE from a
## truncated pipe into a task failure, so every `... | head` is written
## through a temp file. And GLIMPSE2_split_reference appends its own
## _chrom_start_end to whatever output prefix it was handed, so the binary
## reference files carry the chromosome twice: reference_chr22_chr22_1_*.bin.
## Matching on the coordinates alone sidesteps that.

workflow PhaseLigateOnly {
  input {
    String chrom = "chr22"

    File   merged_bcf
    File   merged_csi
    File   chunks_file
    File   ref_bins_tar

    Float  maf_out = 0.001

    String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    String glimpse_docker  = "simrub/glimpse:v2.0.0-27-g0919952_20221207"
  }

  call ReadChunks {
    input:
      chunks = chunks_file,
      docker = bcftools_docker
  }

  scatter (ck in ReadChunks.rows) {
    call Phase {
      input:
        chrom        = chrom,
        chunk_id     = ck[0],
        input_region = ck[2],
        merged_bcf   = merged_bcf,
        merged_csi   = merged_csi,
        ref_bins_tar = ref_bins_tar,
        docker       = glimpse_docker
    }
  }

  call Ligate {
    input:
      chrom      = chrom,
      chunk_bcfs = Phase.phased_bcf,
      chunk_csis = Phase.phased_csi,
      maf_out    = maf_out,
      docker     = glimpse_docker
  }

  output {
    File out_imputed     = Ligate.imputed_bcf
    File out_imputed_csi = Ligate.imputed_csi
    File out_common      = Ligate.common_bcf
    File out_common_csi  = Ligate.common_csi
    Int  out_n_variants  = Ligate.n_variants
    Int  out_n_common    = Ligate.n_common
  }
}


task ReadChunks {
  input {
    File   chunks
    String docker
  }

  command <<<
    set -euo pipefail
    # id, chrom, buffered region, output region; the cM/bp/count columns
    # after those are not needed here.
    awk -F'\t' 'NF>=4 {print $1"\t"$2"\t"$3"\t"$4}' ~{chunks} > rows.tsv
    echo "chunks: $(wc -l < rows.tsv)"
    cat rows.tsv
  >>>

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }

  output {
    Array[Array[String]] rows = read_tsv("rows.tsv")
  }
}


task Phase {
  input {
    String chrom
    String chunk_id
    String input_region
    File   merged_bcf
    File   merged_csi
    File   ref_bins_tar
    String docker
  }

  Int disk_gb = ceil(size(merged_bcf, "GB") * 2 + size(ref_bins_tar, "GB") * 3 + 40)

  command <<<
    set -euo pipefail

    ls -1 /usr/local/bin /usr/bin /opt/*/bin 2>/dev/null > allbin.txt || true
    grep -i glimpse allbin.txt || true
    PHASE=$(command -v GLIMPSE2_phase_static || command -v GLIMPSE2_phase)
    echo "phase=$PHASE"

    ln -s ~{merged_bcf} ./merged.bcf
    ln -s ~{merged_csi} ./merged.bcf.csi

    tar xzf ~{ref_bins_tar}
    # PrepareReference tarred these under refbins/; flatten so the lookup
    # below does not depend on that layout.
    find . -name "reference_*.bin" -exec mv {} . \; 2>/dev/null || true
    ls reference_*.bin > binlist.txt
    echo "reference bins: $(wc -l < binlist.txt)"
    head -3 binlist.txt

    # Match on the buffered coordinates rather than reconstructing the whole
    # filename: split_reference wrote reference_<prefix>_<chrom>_<start>_<end>
    # and the prefix already contained the chromosome.
    IRG_START=$(echo "~{input_region}" | sed 's/.*://; s/-.*//')
    IRG_END=$(echo "~{input_region}" | sed 's/.*-//')
    grep "_${IRG_START}_${IRG_END}\.bin$" binlist.txt > match.txt || true
    BIN=$(head -1 match.txt)
    if [ -z "$BIN" ]; then
      echo "no match for ${IRG_START}_${IRG_END}; available:"
      cat binlist.txt
      exit 1
    fi
    echo "using $BIN"

    "$PHASE" \
      --input-gl ./merged.bcf \
      --reference "$BIN" \
      --output phased_~{chrom}_~{chunk_id}.bcf \
      --threads 4

    bcftools index -f phased_~{chrom}_~{chunk_id}.bcf
    echo "phased ~{chrom} chunk ~{chunk_id}: $(bcftools index -n phased_~{chrom}_~{chunk_id}.bcf) variants"
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "32 GB"
    disks: "local-disk ~{disk_gb} SSD"
    # Phasing thousands of samples over a multi-megabase chunk runs long
    # enough that a reclaim costs more than the discount saves.
    preemptible: 0
    maxRetries: 1
  }

  output {
    File phased_bcf = "phased_~{chrom}_~{chunk_id}.bcf"
    File phased_csi = "phased_~{chrom}_~{chunk_id}.bcf.csi"
  }
}


task Ligate {
  input {
    String      chrom
    Array[File] chunk_bcfs
    Array[File] chunk_csis
    Float       maf_out
    String      docker
  }

  Int disk_gb = ceil(size(chunk_bcfs, "GB") * 4 + 60)

  command <<<
    set -euo pipefail

    LIGATE=$(command -v GLIMPSE2_ligate_static || command -v GLIMPSE2_ligate)
    echo "ligate=$LIGATE"

    mkdir -p work
    while read -r B; do ln -sf "$B" work/$(basename "$B"); done < ~{write_lines(chunk_bcfs)}
    while read -r C; do ln -sf "$C" work/$(basename "$C"); done < ~{write_lines(chunk_csis)}

    # Ligate stitches chunks in coordinate order, and the chunk index in each
    # filename is what establishes that order.
    ls work/phased_~{chrom}_*.bcf | sort -t_ -k3 -n > chunk_list.txt
    echo "chunks to ligate: $(wc -l < chunk_list.txt)"
    cat chunk_list.txt

    "$LIGATE" --input chunk_list.txt --output ~{chrom}_imputed.bcf --threads 4
    bcftools index -f ~{chrom}_imputed.bcf

    NV=$(bcftools index -n ~{chrom}_imputed.bcf)
    echo "$NV" > n_variants.txt
    bcftools query -l ~{chrom}_imputed.bcf > samples.txt
    echo "ligated ~{chrom}: $NV variants, $(wc -l < samples.txt) samples"
    head -3 samples.txt

    # A MAF-filtered copy for downstream PRS work. INFO score is not usable
    # as a filter here: in a BGE design the same reads feed both the
    # likelihoods and the imputation, so it saturates for nearly every site.
    bcftools view -q ~{maf_out}:minor -Ob -o ~{chrom}_common.bcf ~{chrom}_imputed.bcf
    bcftools index -f ~{chrom}_common.bcf

    NC=$(bcftools index -n ~{chrom}_common.bcf)
    echo "$NC" > n_common.txt
    echo "after MAF >= ~{maf_out}: $NC variants"
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: 0
    maxRetries: 1
  }

  output {
    File imputed_bcf = "~{chrom}_imputed.bcf"
    File imputed_csi = "~{chrom}_imputed.bcf.csi"
    File common_bcf  = "~{chrom}_common.bcf"
    File common_csi  = "~{chrom}_common.bcf.csi"
    Int  n_variants  = read_int("n_variants.txt")
    Int  n_common    = read_int("n_common.txt")
  }
}
