process sapin {
    label "sapin"
    publishDir("${params.outdir}/sapin", mode: "copy")

    input:
    tuple(val(sample), path(bam))
    path(reference)

    output:
    path("restable.tsv")

    script:
    """
    sapin -a ${bam} -f ${reference} > restable.tsv
    """
}
