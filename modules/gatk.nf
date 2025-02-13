process gatk_markduplicates {
    label 'gatk'
    tag "$sample_id"
    publishDir "${params.outdir}/gatk_markduplicates", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path ("*_marked_duplicates.bam"), emit: tuple_sample_dedupbam
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
