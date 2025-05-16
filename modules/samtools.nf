process samtools_index {
    label "samtools"
    tag "${meta.id}"
    publishDir("${params.outdir}/bam_indices", mode:"copy")
    
    input:
    tuple(val(meta), path(bam))

    output:
    tuple(val(meta), path(bam), path("*.bai"), emit: tuple_sample_bam_bamindex)

    script:
    """
    samtools index -b ${bam}
    """
}

process samtools_fasta_index {
    label "samtools"
    tag "${fasta}"
    publishDir("${params.outdir}/fasta_indices", mode:"copy")
    
    input:
        path(fasta)

    output:
     tuple(path(fasta), path("${fasta}.fai"), emit: tuple_fasta_fastaindex)

    script:
        """
        samtools faidx ${fasta}
        """
}
