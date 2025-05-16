process sapin {
    label "sapin"
    tag "${meta.id}"
    publishDir("${params.outdir}/sapin", mode: "copy")

    input:
        tuple(val(meta), path(bam))
        path(reference)

    output:
        path("restable.tsv")

    script:
        """
        sapin -a ${bam} -f ${reference} > restable.tsv
        """
}
