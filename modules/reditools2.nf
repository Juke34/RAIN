process reditools2 {
    label "reditools2"
    publishDir("${params.outdir}/reditools", mode: "copy")
    tag "${meta.id}"

    input:
        tuple(val(meta), path(bam), path(bamindex))
        path genome
        val region

    output:
        tuple(val(meta), path("edit_table.txt"), emit: tuple_sample_serial_table)

    script:

    // Set the strand orientation parameter from the library type parameter
    // Terms explained in https://salmon.readthedocs.io/en/latest/library_type.html
    if (meta.libtype in ["ISR", "SR"]) {
        // First-strand oriented
        strand_orientation = "2"
    } else if (meta.libtype in ["ISF", "SF"]) {
        // Second-strand oriented
        strand_orientation = "1"
    } else if (meta.libtype in ["IU", "U"]) {
        // Unstranded
        strand_orientation = "0"
    } else {
        // Unsupported: Pass the library type string so that it's reported in
        // the reditools error message
        strand_orientation = meta.libtype
    }
    
    """
    reditools.py -f ${bam} -r ${genome} -s ${strand_orientation} -o edit_table.txt
    """
}
