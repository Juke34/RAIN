process reditools2 {
    label "reditools2"
    publishDir("${params.outdir}/reditools", mode: "copy")

    input:
    tuple(val(sample), path(bam), path(bamindex))
    path genome
    val region

    output:
    tuple(val(sample), path("edit_table.txt"), emit: tuple_sample_serial_table)

    script:
    """
    reditools.py -f ${bam} -r ${genome} -o edit_table.txt
    """
}
