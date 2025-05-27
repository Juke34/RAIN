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

/*
http://www.htslib.org/doc/samtools-sort.html
Sort alignments by leftmost coordinates 
*/
process samtools_sort_bam {
    label "samtools"
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path ("*_sorted.bam"), emit: tuple_sample_sortedbam

    script:

        """
            if [[ \$(samtools view -H ${bam}| awk '/^@HD/ { for(i=1;i<=NF;i++) if(\$i ~ /^SO:/) print \$i }') == *"coordinate"* ]]; then
                echo "Already sorted by coordinate"
                ln -s ${bam} ${bam.baseName}_sorted.bam
            else
                echo "Not sorted by coordinate"
                samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam  
            fi
        """
}