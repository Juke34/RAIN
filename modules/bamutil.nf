process bamutil_clipoverlap {
    label 'bamutil'
    tag "$sample_id"
    publishDir "${params.outdir}/bamutil_clipoverlap", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path ("*clipoverlap.bam"), emit: tuple_sample_clipoverbam
    path ("*_bamutil_clipoverlap.log"), emit: log

    script:
    """
    bam clipOverlap --storeOrig CG --poolSize 50000000 --in ${bam} --out ${bam}_clipoverlap.bam --stats > ${sample_id}_bamutil_clipoverlap.log
    """

}
