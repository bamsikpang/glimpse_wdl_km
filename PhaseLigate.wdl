version 1.0

## PhaseLigate.wdl
##
## Merges the per-sample genotype likelihoods into one BCF per chromosome,
## runs GLIMPSE2_phase over each chunk, and ligates the chunks back into a
## whole-chromosome imputed callset.
##
## Chunk boundaries come from the chunks_chrN.txt that PrepareReference wrote,
## so the scatter width is whatever GLIMPSE2_chunk decided: 433 chunks across
## chr1-chr21, from 7 on chr21 to 36 on chr1.
##
## Four things learned from the chr22 run are baked in here.
##
## The bcftools image ships BusyBox coreutils, so `split -d` is unavailable
## and the batching is done with awk. `set -o pipefail` turns a SIGPIPE from
## a truncated pipe into a task failure, so every `... | head` writes through
## a temp file first. GLIMPSE2_split_reference appends its own
## _chrom_start_end to whatever output prefix it received, leaving filenames
## that carry the chromosome twice, so the binary reference lookup matches on
## coordinates alone. And the CRAMs carry their original delivery IDs while
## the phenotype tables use a normalized form, so samples are renamed after
## the merge.

workflow PhaseLigate {
  input {
    Array[String] chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                            "chr8","chr9","chr10","chr11","chr12","chr13",
                            "chr14","chr15","chr16","chr17","chr18","chr19",
                            "chr20","chr21"]

    String bcflist_dir              # holds chrN.txt and chrN.csi.txt
    String panel_dir                # holds chunks_chrN.txt and refbins_chrN.tar.gz

    Int    merge_batch = 400
    Float  maf_out     = 0.001       # post-imputation MAF floor for the filtered copy
    Float  info_min    = 0.8        # IMPUTE-style quality floor

    String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
    String glimpse_docker  = "simrub/glimpse:v2.0.0-27-g0919952_20221207"
  }

  scatter (c in chroms) {

    call MergeGL {
      input:
        chrom      = c,
        bcf_list   = bcflist_dir + "/" + c + ".txt",
        csi_list   = bcflist_dir + "/" + c + ".csi.txt",
        batch_size = merge_batch,
        docker     = bcftools_docker
    }

    call ReadChunks {
      input:
        chunks = panel_dir + "/chunks_" + c + ".txt",
        docker = bcftools_docker
    }

    scatter (ck in ReadChunks.rows) {
      call Phase {
        input:
          chrom        = c,
          chunk_id     = ck[0],
          input_region = ck[2],
          merged_bcf   = MergeGL.merged_bcf,
          merged_csi   = MergeGL.merged_csi,
          ref_bins_tar = panel_dir + "/refbins_" + c + ".tar.gz",
          docker       = glimpse_docker
      }
    }

    call Ligate {
      input:
        chrom      = c,
        chunk_bcfs = Phase.phased_bcf,
        chunk_csis = Phase.phased_csi,
        maf_out    = maf_out,
        info_min   = info_min,
        docker     = glimpse_docker
    }
  }

  output {
    Array[File] out_imputed     = Ligate.imputed_bcf
    Array[File] out_imputed_csi = Ligate.imputed_csi
    Array[File] out_common      = Ligate.common_bcf
    Array[File] out_common_csi  = Ligate.common_csi
    Array[Int]  out_n_variants  = Ligate.n_variants
    Array[Int]  out_n_common    = Ligate.n_common
    Array[File] out_merge_log   = MergeGL.merge_log
    Array[File] out_rename_map  = MergeGL.rename_map
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


task MergeGL {
  input {
    String chrom
    File   bcf_list
    File   csi_list
    Int    batch_size
    String docker
  }

  Array[File] bcfs = read_lines(bcf_list)
  Array[File] csis = read_lines(csi_list)

  Int disk_gb = ceil(size(bcfs, "GB") * 3 + 60)

  command <<<
    set -euo pipefail

    cat ~{write_lines(bcfs)} > bcf_paths.txt
    cat ~{write_lines(csis)} > csi_paths.txt
    echo "bcfs: $(wc -l < bcf_paths.txt)  csis: $(wc -l < csi_paths.txt)"

    # An index has to sit beside its data file, and Cromwell localizes each
    # input into its own directory. Rebuild a flat working directory of
    # symlinks so htslib can pair them.
    mkdir -p work
    while read -r B; do ln -sf "$B" work/$(basename "$B"); done < bcf_paths.txt
    while read -r C; do ln -sf "$C" work/$(basename "$C"); done < csi_paths.txt
    ls work/*.bcf | sort > work_bcfs.txt
    echo "linked: $(wc -l < work_bcfs.txt)"

    # bcftools merge holds every input open at once and 2,420 exceeds the
    # usual 1024-descriptor ceiling, so merge in batches and then merge the
    # batches. BusyBox split has no -d, hence awk.
    mkdir -p batches
    awk -v n=~{batch_size} '{
      printf("%s\n", $0) > sprintf("batches/part_%03d", int((NR-1)/n))
    }' work_bcfs.txt
    N_BATCH=$(ls batches/part_* | wc -l)
    echo "batches: $N_BATCH"

    for P in batches/part_*; do
      OUT=batches/$(basename $P).bcf
      bcftools merge -m none -l "$P" -Ob -o "$OUT"
      bcftools index -f "$OUT"
      echo "  $(basename $P): $(bcftools query -l $OUT | wc -l) samples"
    done

    ls batches/part_*.bcf | sort > batch_list.txt
    if [ "$N_BATCH" -eq 1 ]; then
      cp $(head -1 batch_list.txt) merged_raw.bcf
    else
      bcftools merge -m none -l batch_list.txt -Ob -o merged_raw.bcf
    fi
    bcftools index -f merged_raw.bcf

    # The CRAMs carried their original delivery IDs (ABN-K1-HC-00001) while
    # the manifest and phenotype tables use the normalized form sitting in
    # the filename (ABN-KR1-HC-00001). bcftools merge preserves input order,
    # so the sorted filename list gives the replacement names directly.
    awk -F/ '{print $NF}' work_bcfs.txt | cut -d. -f1 > newnames.txt
    bcftools query -l merged_raw.bcf > oldnames.txt
    paste oldnames.txt newnames.txt > rename_~{chrom}.map
    echo "rename map (first 3):"
    head -3 rename_~{chrom}.map
    echo "old $(wc -l < oldnames.txt) / new $(wc -l < newnames.txt)"

    bcftools reheader -s newnames.txt -o merged_~{chrom}.bcf merged_raw.bcf
    bcftools index -f merged_~{chrom}.bcf

    bcftools query -l merged_~{chrom}.bcf > samplenames.txt
    NS=$(wc -l < samplenames.txt)
    NV=$(bcftools index -n merged_~{chrom}.bcf)
    echo "merged ~{chrom}: $NS samples, $NV variants" | tee merge_~{chrom}.log
    head -3 samplenames.txt

    # GLIMPSE2 reads PL or GL; confirm the field survived the merge.
    bcftools view -h merged_~{chrom}.bcf > hdr.txt
    grep -E "FORMAT=<ID=(PL|GL)" hdr.txt || echo "WARNING: no PL/GL in FORMAT"
  >>>

  runtime {
    docker: docker
    cpu: 4
    memory: "16 GB"
    disks: "local-disk ~{disk_gb} SSD"
    # Localizing 4,840 per-sample files ran 2h40m for chr22 alone; a reclaim
    # partway through would throw all of that away.
    preemptible: 0
    maxRetries: 1
  }

  output {
    File merged_bcf = "merged_~{chrom}.bcf"
    File merged_csi = "merged_~{chrom}.bcf.csi"
    File merge_log  = "merge_~{chrom}.log"
    File rename_map = "rename_~{chrom}.map"
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

    # Match on the buffered coordinates rather than rebuilding the filename:
    # split_reference wrote reference_<prefix>_<chrom>_<start>_<end> and the
    # prefix already ended in the chromosome, so the name has it twice.
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
    # chr22 chunks took 4.8 to 7.3 hours at 2,420 samples, far past the point
    # where preemption pays for itself.
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
    Float       info_min
    String      docker
  }

  Float maf_max = 1.0 - maf_out
  Int   disk_gb = ceil(size(chunk_bcfs, "GB") * 4 + 60)

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

    # GLIMPSE2 derives INFO/AF from the dosages and INFO/INFO as an
    # IMPUTE-style quality score. Filtering on GT instead would throw away
    # low-frequency sites whose hard calls all collapse to 0/0. The quality
    # score was saturated in a 516-sample subset but separates properly
    # across all 2,420 (chr22 mean 0.861), so it earns a place here.
    bcftools view \
      -i 'INFO/AF >= ~{maf_out} && INFO/AF <= ~{maf_max} && INFO/INFO >= ~{info_min}' \
      -Ob -o ~{chrom}_common.bcf ~{chrom}_imputed.bcf
    bcftools index -f ~{chrom}_common.bcf

    NC=$(bcftools index -n ~{chrom}_common.bcf)
    echo "$NC" > n_common.txt
    echo "after MAF >= ~{maf_out} and INFO >= ~{info_min}: $NC variants"
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
