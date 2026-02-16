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
                ln -s \$(realpath ${bam}) ${bam.baseName}_sorted.bam
            else
                echo "Not sorted by coordinate"
                samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam  
            fi
        """
}

/*
 * Split BAM file into mapped and unmapped reads
 */
process samtools_split_mapped_unmapped {
    label "samtools"
    tag "${meta.id}"
    publishDir("${params.outdir}/split_bam", mode:"copy", pattern: "*_{mapped,unmapped}.bam")
    
    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*_mapped.bam"), emit: mapped_bam
        tuple val(meta), path("*_unmapped.bam"), emit: unmapped_bam

    script:
        """
        # Extract mapped reads (SAM flag -F 4: exclude unmapped)
        samtools view -@ ${task.cpus} -b -F 4 ${bam} > ${bam.baseName}_mapped.bam
        
        # Extract unmapped reads (SAM flag -f 4: include only unmapped)
        samtools view -@ ${task.cpus} -b -f 4 ${bam} > ${bam.baseName}_unmapped.bam
        """
}

/*
 * Convert BAM to FASTQ format
 * Extracts reads from BAM and converts to FASTQ for re-alignment
 */
process convert_to_fastq {
    label "samtools"
    tag "${meta.id}"
    publishDir("${output_dir}", mode:"copy", pattern: "*.fastq*")
    
    input:
        tuple val(meta), path(bam)
        val output_dir

    output:
        tuple val(meta), path("*.fastq.gz"), emit: fastq_reads

    script:
        def output_base = "${bam.baseName}"
        """
        # Check if BAM contains paired-end or single-end reads
        if samtools view -c -f 1 ${bam} | grep -q "^0\$"; then
            # Single-end reads
            samtools fastq -@ ${task.cpus} ${bam} | gzip > ${output_base}.fastq.gz
        else
            # Paired-end reads
            samtools fastq -@ ${task.cpus} -1 ${output_base}_R1.fastq.gz -2 ${output_base}_R2.fastq.gz -s ${output_base}_single.fastq.gz ${bam}
        fi
        """
}

/*
 * Merge two BAM files
 * Combines aligned reads from two BAM files into a single output BAM
 */
process samtools_merge_bams {
    label "samtools"
    tag "${meta.id}"
    publishDir("${params.outdir}/merged_bam", mode:"copy")
    
    input:
        tuple val(meta), path(bam1), path(bam2)

    output:
        tuple val(meta), path("*_merged.bam"), emit: merged_bam

    script:
        """
        samtools merge -@ ${task.cpus} ${meta.id}_merged.bam ${bam1} ${bam2}
        """
}