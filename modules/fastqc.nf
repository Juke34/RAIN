process fastqc {
    label 'fastqc'
    tag "${meta.uid}"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        val (suffix)

    output:
        path ("fastqc_${meta.uid}_logs_${suffix}")

    script:
        """
        mkdir fastqc_${meta.uid}_logs_${suffix}
        fastqc -t ${task.cpus} -o fastqc_${meta.uid}_logs_${suffix} -q ${reads}
        """

}
