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
    if (params.library_type in ["ISR", "SR"]) {
        // First-strand oriented
        strand_orientation = "2"
    } else if (params.library_type in ["ISF", "SF"]) {
        // Second-strand oriented
        strand_orientation = "1"
    } else if (params.library_type in ["IU", "U"]) {
        // Unstranded
        strand_orientation = "0"
    } else {
        // Unsupported: Pass the library type string so that it's reported in
        // the reditools error message
        strand_orientation = params.library_type
    }
    """
    reditools.py -f ${bam} -r ${genome} -s ${strand_orientation} -o edit_table.txt
    """
}
