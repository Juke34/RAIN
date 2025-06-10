process jacusa2 {
    label "jacusa2"
    tag "${meta.id}"
    publishDir("${params.outdir}/jacusa2", mode: "copy")

    input:
        tuple(val(meta), path(bam), path(bamindex))
        tuple(path(genome), path(fastaindex))

    output:
        tuple(val(meta), path("*.site_edits_jacusa2.txt"), emit: tuple_sample_jacusa2_table)
        path("*.filtered")

    script:
        base_name = bam.BaseName
        // Set the strand orientation parameter from the library type parameter
        // Terms explained in https://salmon.readthedocs.io/en/latest/library_type.html
        if (meta.libtype in ["ISR", "SR"]) {
            // First-strand oriented
            jacusa_strand_param = "FR_SECONDSTRAND"
        } else if (meta.libtype in ["ISF", "SF"]) {
            // Second-strand oriented
            jacusa_strand_param = "RF_FIRSTSTRAND"
        } else if (meta.libtype in ["IU", "U"]) {
            // Unstranded
            jacusa_strand_param = "UNSTRANDED"
        } else {
            // Unsupported: Pass the library type string so that it's reported in
            // the Jacusa2 error message
            jacusa_strand_param = meta.libtype
        }
        """
        java -Xmx${task.memory.toMega()}M -jar /usr/local/bin/JACUSA_v2.0.4.jar \
            call-1 \
            -A \
            -f V \
            -p ${task.cpus} \
            -r ${base_name}.site_edits_jacusa2.txt \
            -c 1 \
            -s \
            -R ${genome} \
            -P ${jacusa_strand_param} \
            ${bam}                                      
        """
}