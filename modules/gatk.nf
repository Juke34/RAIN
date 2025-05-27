process gatk_markduplicates {
    label 'gatk'
    tag "${meta.id}"
    publishDir "${params.outdir}/gatk_markduplicates", mode: 'copy'

    input:
      tuple val(meta), path(bam)

    output:
      tuple val(meta), path ("*_marked_duplicates.bam"), emit: tuple_sample_dedupbam
      path ("*_marked_dup_metrics.txt") , emit: log

    script:

      """
      gatk MarkDuplicates \
        -I ${bam} \
        -O ${bam.baseName}_marked_duplicates.bam \
        -M ${bam.baseName}_marked_dup_metrics.txt \
        --REMOVE_DUPLICATES
      """
}
