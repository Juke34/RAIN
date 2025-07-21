process genome_recoder {
    label "genome_recoder"
    publishDir("${params.outdir}/recoded_genome", mode: "copy")
    tag "$genome"

    input:
        path(genome)

    output:
        path("*.recoded_AG.fasta")

    script:
        base_name = genome.BaseName
        """
        genome_recoder.py --genome ${genome} --base A --replacement G --prefix ${base_name}.recoded_
        """
}