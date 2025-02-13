process fastqc {
    label 'fastqc'
    tag "$sample_id"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val (suffix)

    output:
    path ("fastqc_${sample_id}_logs_${suffix}")

    script:
    """
    mkdir fastqc_${sample_id}_logs_${suffix}
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs_${suffix} -q ${reads}
    """

}
