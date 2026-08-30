version 1.0

## ComputeGL.wdl
##
## Computes genotype likelihoods at reference-panel sites for each sample,
## feeding GLIMPSE2_phase.
##
## Scatter is over samples, not over sample-by-chromosome pairs. A CRAM is
## 2.78 GB and the earlier hand-run pilot showed localization dominating
## runtime; splitting by chromosome would move 148 TB instead of 6.7 TB while
## the per-chromosome pileups themselves are cheap once the file is local.
##
## The 88 site files (22 chromosomes x 4) arrive as one 198 MB tarball. The
## bcftools image carries no gsutil, so fetching them inside the task is not
## available, and 88 separate File inputs would make chromosome-to-file
## pairing fragile.
##
## Default chroms is chr22 alone. That chromosome already has a hand-run
## reference point (873 s per sample, 222,313 sites) to check this pipeline
## against, and it costs 15 minutes per sample instead of hours. Override
## chroms in the inputs JSON to run the full set.

workflow ComputeGL {
  input {
    File   sample_manifest          # sid, cram, crai (tab-separated, with header)
    File   sites_tar                # sites_chrN.{vcf.gz,vcf.gz.tbi,tsv.gz,tsv.gz.tbi}
    File   ref_fasta
    File   ref_fasta_index

    Array[String] chroms = ["chr22"]

    Int    mapq_min  = 30
    Int    baseq_min = 20

    String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
  }

  call ReadManifest {
    input:
      manifest = sample_manifest,
      docker   = bcftools_docker
  }

  scatter (row in ReadManifest.rows) {
    call SampleGL {
      input:
        sid             = row[0],
        cram            = row[1],
        crai            = row[2],
        sites_tar       = sites_tar,
        chroms          = chroms,
        ref_fasta       = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        mapq_min        = mapq_min,
        baseq_min       = baseq_min,
        docker          = bcftools_docker
    }
  }

  output {
    Array[File] out_gl_tars = SampleGL.gl_tar
    Array[File] out_timings = SampleGL.timing
  }
}


task ReadManifest {
  input {
    File   manifest
    String docker
  }

  command <<<
    set -euo pipefail
    tail -n +2 ~{manifest} | awk -F'\t' 'NF>=3 {print $1"\t"$2"\t"$3}' > rows.tsv
    echo "samples: $(wc -l < rows.tsv)"
    head -2 rows.tsv
  >>>

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    # Seconds of work; preemption has no window in which to bite.
    preemptible: 3
  }

  output {
    Array[Array[String]] rows = read_tsv("rows.tsv")
  }
}


task SampleGL {
  input {
    String sid
    File   cram
    File   crai
    File   sites_tar
    Array[String] chroms
    File   ref_fasta
    File   ref_fasta_index
    Int    mapq_min
    Int    baseq_min
    String docker
  }

  # CRAM, reference, unpacked sites, per-chromosome BCFs, output tarball.
  Int disk_gb = ceil(size(cram, "GB") + size(ref_fasta, "GB") + size(sites_tar, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    # A bare "command not found" surfaces as exit 127 with no indication of
    # which tool was absent, so name them up front.
    for t in bcftools bgzip tabix; do
      command -v $t || echo "MISSING: $t"
    done

    # htslib resolves CRAM reference sequences through REF_CACHE, and when
    # that is unset it will try the EBI over the network one checksum at a
    # time. Point it at a local directory so -f supplies the sequence and
    # repeated chromosomes reuse what was already read.
    mkdir -p ref_cache
    export REF_CACHE=$PWD/ref_cache/%2s/%2s/%s
    export REF_PATH=$PWD/ref_cache/%2s/%2s/%s

    tar xzf ~{sites_tar}
    echo "site files unpacked: $(ls sites_chr*.vcf.gz | wc -l)"

    # htslib expects an index beside its data file.
    ln -s ~{cram} ./sample.cram
    ln -s ~{crai} ./sample.cram.crai

    mkdir -p gl
    : > timing.txt

    for CHR in ~{sep=' ' chroms}; do
      T0=$(date +%s)

      # -I skips indels, -E recomputes BAQ, -T restricts the pileup to panel
      # sites. call -C alleles forces the panel's REF/ALT rather than letting
      # bcftools infer them from the reads, which is what GLIMPSE2 assumes.
      bcftools mpileup \
        -f ~{ref_fasta} \
        -I -E \
        -a 'FORMAT/DP' \
        -q ~{mapq_min} -Q ~{baseq_min} \
        -T ./sites_${CHR}.vcf.gz \
        -r ${CHR} \
        ./sample.cram -Ou \
        | bcftools call -Aim -C alleles -T ./sites_${CHR}.tsv.gz \
            -Ob -o gl/~{sid}.${CHR}.bcf
      bcftools index -f gl/~{sid}.${CHR}.bcf

      T1=$(date +%s)
      N_OUT=$(bcftools index -n gl/~{sid}.${CHR}.bcf)
      N_SITE=$(bcftools index -n ./sites_${CHR}.vcf.gz)
      echo -e "~{sid}\t${CHR}\t$((T1-T0))\t${N_OUT}\t${N_SITE}" >> timing.txt
      echo "~{sid} ${CHR}: $((T1-T0))s, ${N_OUT}/${N_SITE} sites called"
    done

    tar czf ~{sid}.gl.tar.gz -C gl .

    echo "=== summary ==="
    ls -lh ~{sid}.gl.tar.gz
    cat timing.txt
    awk -F'\t' '{s+=$3; o+=$4; p+=$5}
      END {printf "total %ds | called %d / panel %d (%.1f%%)\n", s, o, p, 100*o/p}' timing.txt
  >>>

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GB"
    disks: "local-disk ~{disk_gb} SSD"
    # Whole-genome pileup is long enough that preemption costs more in
    # wasted restarts than it saves; PrepareReference lost three hours to
    # two reclaims on a 1.5 h task.
    preemptible: 0
    maxRetries: 1
  }

  output {
    File gl_tar = "~{sid}.gl.tar.gz"
    File timing = "timing.txt"
  }
}
