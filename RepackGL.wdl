version 1.0

## RepackGL.wdl
##
## Turns per-sample GL tarballs into per-chromosome BCFs.
##
## ComputeGL emitted one tarball per sample holding every chromosome. Merging
## chromosome by chromosome from that layout would pull all 482 GB once per
## chromosome, 10.6 TB in total, to extract a few megabytes each time. Unpacking
## once here costs a single 482 GB pass and leaves MergeGL fetching ~22 GB per
## chromosome.
##
## Samples were run in two batches (chr1-21 and chr22 separately), so each
## sample has two tarballs. Both are unpacked into the same directory.

workflow RepackGL {
  input {
    File   gl_tar_list_main     # one gs:// URL per line, chr1-21 tarballs
    File   gl_tar_list_chr22    # one gs:// URL per line, chr22 tarballs
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
  }

  Array[File] tars_main  = read_lines(gl_tar_list_main)
  Array[File] tars_chr22 = read_lines(gl_tar_list_chr22)

  scatter (i in range(length(tars_main))) {
    call Unpack {
      input:
        tar_main  = tars_main[i],
        tar_chr22 = tars_chr22[i],
        docker    = bcftools_docker
    }
  }

  output {
    Array[Array[File]] out_bcfs = Unpack.bcfs
    Array[Array[File]] out_csis = Unpack.csis
    Array[String]      out_sids = Unpack.sid
  }
}


task Unpack {
  input {
    File   tar_main
    File   tar_chr22
    String docker
  }

  Int disk_gb = ceil((size(tar_main, "GB") + size(tar_chr22, "GB")) * 4 + 10)

  command <<<
    set -euo pipefail

    tar xzf ~{tar_main}
    tar xzf ~{tar_chr22}

    # ComputeGL wrote these as <sid>.<chrom>.bcf, so the sample name is
    # whatever precedes the first dot on any of them.
    SID=$(ls *.bcf | head -1 | cut -d. -f1)
    echo "$SID" > sid.txt

    echo "sample: $SID"
    echo "bcfs: $(ls *.bcf | wc -l)"
    ls *.bcf
  >>>

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    # Under a minute of work; a reclaim costs nothing to redo.
    preemptible: 3
    maxRetries: 2
  }

  output {
    Array[File] bcfs = glob("*.bcf")
    Array[File] csis = glob("*.bcf.csi")
    String      sid  = read_string("sid.txt")
  }
}
