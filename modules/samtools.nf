

process samtools_sam_to_bam {
    label 'samtools'
    tag "$sample"
    publishDir "${params.outdir}/Hisat2_alignments", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path ("*.bam"), emit: tuple_sample_bam

    script:

    if (params.single_end){
    """
        samtools view -@ ${task.cpus} ${sam} -b -o ${sam.baseName}.bam 
    """
    } else {
    """
        samtools view -@ ${task.cpus} ${sam} -b -o ${sam.baseName}.bam 
    """
    }

}

process samtools_sort {
    label 'samtools'
    tag "$sample"
    publishDir "${params.outdir}/Hisat2_alignments_sorted", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path ("*_sorted.bam"), emit: tuple_sample_sortedbam

    script:

    if (params.single_end){
    """
        samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam 
    """
    } else {
    """
        samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam  
    """
    }

}

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
