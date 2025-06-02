process sapin {
    label "sapin"
    tag "${meta.id}"
    publishDir("${params.outdir}/sapin", mode: "copy")

    input:
        tuple(val(meta), path(bam))
        path(reference)

    output:
        path("*.site_edits_sapin.tsv")

    script:
        base_name = bam.BaseName
        """
        sapin -a ${bam} -f ${reference} > ${base_name}.site_edits_sapin.tsv
        """
}
