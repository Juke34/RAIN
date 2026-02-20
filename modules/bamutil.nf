process bamutil_clipoverlap {
    label 'bamutil'
    tag "${meta.uid}"
    publishDir "${params.outdir}/bamutil_clipoverlap", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path ("*clipoverlap.bam"), emit: tuple_sample_clipoverbam
    path ("*_bamutil_clipoverlap.log"), emit: log

    script:

        """
        bam clipOverlap --storeOrig CG --poolSize 5000000 --in ${bam} --out ${bam}_clipoverlap.bam --stats --excludeFlags 0x37 > ${meta.uid}_bamutil_clipoverlap.log
        """

}
