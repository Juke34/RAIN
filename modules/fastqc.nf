process fastqc {
    label 'fastqc'
    tag "${meta.id}"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        val (suffix)

    output:
        path ("fastqc_${meta.id}_logs_${suffix}")

    script:
        """
        mkdir fastqc_${meta.id}_logs_${suffix}
        fastqc -t ${task.cpus} -o fastqc_${meta.id}_logs_${suffix} -q ${reads}
        """

}
