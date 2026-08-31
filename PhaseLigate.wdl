version 1.0

## PhaseLigate.wdl
##
## Third and last Cromwell stage: merge the per-sample genotype likelihoods
## into one BCF per chromosome, run GLIMPSE2_phase over each chunk, and
## ligate the chunks back into a whole-chromosome imputed callset.
##
## Chunk boundaries come from the chunks_chrN.txt that PrepareReference wrote,
## so the scatter width is whatever GLIMPSE2_chunk decided: 440 chunks across
## the 22 autosomes, from 7 on chr22 to 36 on chr1.

workflow PhaseLigate {
  input {
    Array[String] chroms = ["chr22"]

    String bcflist_dir              # holds chrN.txt and chrN.csi.txt
    String panel_dir                # holds chunks_chrN.txt and refbins_chrN.tar.gz

    Int    merge_batch = 400
    Float  maf_out     = 0.001      # post-imputation MAF floor for the filtered copy

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
          chrom         = c,
          chunk_id      = ck[0],
          input_region  = ck[2],
          merged_bcf    = MergeGL.merged_bcf,
          merged_csi    = MergeGL.merged_csi,
          ref_bins_tar  = panel_dir + "/refbins_" + c + ".tar.gz",
          docker        = glimpse_docker
      }
    }

    call Ligate {
      input:
        chrom        = c,
        chunk_bcfs   = Phase.phased_bcf,
        chunk_csis   = Phase.phased_csi,
        maf_out      = maf_out,
        docker       = glimpse_docker
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

    # An index has to sit beside its data file, and Cromwell scatters each
    # localized input into its own directory. Rebuild a flat working
    # directory of symlinks so htslib can pair them.
    mkdir -p work
    while read -r B; do ln -sf "$B" work/$(basename "$B"); done < bcf_paths.txt
    while read -r C; do ln -sf "$C" work/$(basename "$C"); done < csi_paths.txt
    ls work/*.bcf | sort > work_bcfs.txt
    echo "linked: $(wc -l < work_bcfs.txt)"

    # bcftools merge holds every input open at once; 2,420 exceeds the usual
    # 1024-descriptor ceiling, so merge in batches and then merge the batches.
    mkdir -p batches
    split -l ~{batch_size} -d -a 3 work_bcfs.txt batches/part_
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

    # The CRAMs carried their original delivery IDs (ABN-K11-PT-0001) while
    # the manifest and phenotype tables use the normalized form sitting in
    # the filename (ABN-KR11-PT-00001). bcftools merge preserves input order,
    # so the sorted filename list gives the replacement names directly.
    awk -F/ '{print $NF}' work_bcfs.txt | cut -d. -f1 > newnames.txt
    paste <(bcftools query -l merged_raw.bcf) newnames.txt > rename_~{chrom}.map
    echo "rename map (first 3):"
    head -3 rename_~{chrom}.map
    echo "names to apply: $(wc -l < newnames.txt)"

    bcftools reheader -s newnames.txt -o merged_~{chrom}.bcf merged_raw.bcf
    bcftools index -f merged_~{chrom}.bcf

    NS=$(bcftools query -l merged_~{chrom}.bcf | wc -l)
    NV=$(bcftools index -n merged_~{chrom}.bcf)
    echo "merged ~{chrom}: $NS samples, $NV variants" | tee merge_~{chrom}.log
    bcftools query -l merged_~{chrom}.bcf | head -3

    # GLIMPSE2 reads PL or GL; confirm the field survived the merge.
    bcftools view -h merged_~{chrom}.bcf | grep -E "FORMAT=<ID=(PL|GL)" \
      || echo "WARNING: no PL/GL in FORMAT"
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

    ls -1 /usr/local/bin /usr/bin /opt/*/bin 2>/dev/null | grep -i glimpse || true
    PHASE=$(command -v GLIMPSE2_phase_static || command -v GLIMPSE2_phase)
    echo "phase=$PHASE"

    ln -s ~{merged_bcf} ./merged.bcf
    ln -s ~{merged_csi} ./merged.bcf.csi

    tar xzf ~{ref_bins_tar}
    # PrepareReference tarred these under refbins/; flatten so the name match
    # below does not depend on that layout.
    find . -name "reference_*.bin" -exec mv {} . \; 2>/dev/null || true
    ls reference_*.bin | head -3

    # split_reference named each file after the region it covers, so the
    # buffered coordinates from chunks.txt recover the right one.
    IRG_START=$(echo "~{input_region}" | sed 's/.*://; s/-.*//')
    IRG_END=$(echo "~{input_region}" | sed 's/.*-//')
    BIN=$(ls reference_~{chrom}_${IRG_START}_${IRG_END}.bin 2>/dev/null || true)
    if [ -z "$BIN" ]; then
      echo "no exact match for ${IRG_START}_${IRG_END}; available:"
      ls reference_~{chrom}_*.bin
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
    echo "ligated ~{chrom}: $NV variants, $(bcftools query -l ~{chrom}_imputed.bcf | wc -l) samples"

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
