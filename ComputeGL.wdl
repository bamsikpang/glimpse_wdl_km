version 1.0

## ComputeGL.wdl
##
## Computes genotype likelihoods at the reference panel sites for each sample,
## as input to GLIMPSE2_phase.
##
## The scatter is over samples rather than over sample-by-chromosome pairs.
## A CRAM is 2.78 GB and the pilot showed localization dominating runtime;
## splitting by chromosome would move 148 TB instead of 6.7 TB and gain
## nothing, since the per-chromosome pileups are cheap once the file is on
## local disk.
##
## The 88 site files (22 chromosomes x 4) arrive as one 198 MB tarball. The
## bcftools image has no gsutil, so fetching them inside the task is not an
## option, and passing 88 separate File inputs would make the chromosome-to-
## file pairing fragile.

workflow ComputeGL {
  input {
    File   sample_manifest          # sid, cram, crai (tab-separated, with header)
    File   sites_tar                # sites_chrN.{vcf.gz,vcf.gz.tbi,tsv.gz,tsv.gz.tbi}
    File   ref_fasta
    File   ref_fasta_index

    Array[String] chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                            "chr8","chr9","chr10","chr11","chr12","chr13",
                            "chr14","chr15","chr16","chr17","chr18","chr19",
                            "chr20","chr21","chr22"]

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
    wc -l rows.tsv
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

  # CRAM, reference, unpacked sites, 22 small BCFs, and the output tarball.
  Int disk_gb = ceil(size(cram, "GB") + size(ref_fasta, "GB") + size(sites_tar, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    for t in bcftools bgzip tabix; do
      command -v $t || echo "MISSING: $t"
    done

    # Point htslib at the local FASTA rather than letting it try to fetch
    # reference sequences over the network.
    export REF_PATH=/dev/null
    export REF_CACHE=/dev/null

    tar xzf ~{sites_tar}
    ls sites_chr*.vcf.gz | wc -l

    # htslib expects an index beside its data file.
    ln -s ~{cram} ./sample.cram
    ln -s ~{crai} ./sample.cram.crai

    mkdir -p gl
    : > timing.txt

    for CHR in ~{sep=' ' chroms}; do
      T0=$(date +%s)

      # -I skips indels, -E recomputes BAQ, -T restricts the pileup to panel
      # sites. call -C alleles forces the panel's REF/ALT instead of letting
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
      N=$(bcftools index -n gl/~{sid}.${CHR}.bcf)
      echo -e "~{sid}\t${CHR}\t$((T1-T0))\t${N}" >> timing.txt
      echo "~{sid} ${CHR}: $((T1-T0))s, ${N} sites"
    done

    tar czf ~{sid}.gl.tar.gz -C gl .
    ls -lh ~{sid}.gl.tar.gz
    echo "=== timing ==="
    cat timing.txt
    awk -F'\t' '{s+=$3} END {print "total:", s, "sec"}' timing.txt
  >>>

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GB"
    disks: "local-disk ~{disk_gb} SSD"
    # Whole-genome pileup runs long enough that preemption wastes real work,
    # but a single retry is cheaper than paying full price for every shard.
    preemptible: 1
    maxRetries: 2
  }

  output {
    File gl_tar = "~{sid}.gl.tar.gz"
    File timing = "timing.txt"
  }
}
