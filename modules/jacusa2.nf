process jacusa2 {
    label "jacusa2"
    tag "${meta.id}"
    publishDir("${params.outdir}/jacusa2", mode: "copy")

    input:
        tuple(val(meta), path(bam), path(bamindex))
        tuple(path(genome), path(fastaindex))

    output:
        path("*.site_edits_jacusa2.txt")
        path("*.filtered")
        // tuple(sample, path("filtered_output.txt", emit: tuple_sample_jacusa2_table))

    script:
        base_name = bam.BaseName
        
        """
        java -Xmx${task.memory.toMega()}M -jar /usr/local/bin/JACUSA_v2.0.4.jar call-1 -a D -f V -p ${task.cpus} -r ${base_name}.site_edits_jacusa2.txt -c 1 -s -R ${genome} ${bam}
        """
}