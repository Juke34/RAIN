process samtools_index {
    label "samtools"
    tag "$sample"
    publishDir("${params.outdir}/bam_indices", mode:"copy")
    
    input:
    tuple(val(sample), path(bam))

    output:
    tuple(val(sample), path(bam), path("*.bai"), emit: tuple_sample_bam_bamindex)

    script:
    """
    samtools index -b ${bam}
    """
}

process samtools_fasta_index {
    label "samtools"
    tag "genome"
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
