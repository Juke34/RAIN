process reditools2 {
    label "reditools2"
    publishDir("${params.outdir}/reditools2", mode: "copy")
    tag "${meta.uid}"

    input:
        tuple(val(meta), path(bam), path(bamindex))
        path genome
        val region

    output:
        tuple(val(meta), path("*site_edits_reditools2.txt"), emit: tuple_sample_serial_table)

    script:

    // Set the strand orientation parameter from the library type parameter
    // Terms explained in https://salmon.readthedocs.io/en/latest/library_type.html
    if (meta.strandedness in ["ISR", "SR"]) {
        // First-strand oriented
        strand_orientation = "2"
    } else if (meta.strandedness in ["ISF", "SF"]) {
        // Second-strand oriented
        strand_orientation = "1"
    } else if (meta.strandedness in ["IU", "U"]) {
        // Unstranded
        strand_orientation = "0"
    } else {
        // Unsupported: Pass the library type string so that it's reported in
        // the reditools error message
        print("invalid strand \"${meta.strandedness}\"")
        strand_orientation = 0
    }
    base_name = bam.BaseName
    
    """
    reditools.py -f ${bam} -r ${genome} -s ${strand_orientation} -o ${base_name}.site_edits_reditools2.txt
    """
}
