process jacusa2 {
    label "jacusa2"
    publishDir("${params.outdir}/jacusa2", mode: "copy")

    input:
    tuple(val(sample), path(bam), path(bamindex))
    tuple(path(genome), path(fastaindex))

    output:
    path("*.txt")
    path("*.filtered")
    // tuple(sample, path("filtered_output.txt", emit: tuple_sample_jacusa2_table))

    script:
    """
    java -Xmx${task.memory.toMega()}M -jar /usr/local/bin/JACUSA_v2.0.4.jar call-1 -a D -f V -p ${task.cpus} -r jacusa_out.txt -c 1 -s -R ${genome} ${bam}
    """
}