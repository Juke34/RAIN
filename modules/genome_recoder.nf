process genome_recoder {
    label "genome_recoder"
    publishDir("${params.outdir}/recoded_genome", mode: "copy")
    tag "${meta.id}"

    input:
        tuple(val(meta), path(genome_path))

    output:
        tuple(val(meta), path("*.recoded_AG.fasta"))

    script:
        base_name = genome_path.BaseName
        """
        genome_recoder.py
            --genome ${genome_path} \
            --base A \
            --replacement G \
            --prefix ${base_name}.recoded_
        """
}