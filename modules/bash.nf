/*
Here are described all processes related to bash
*/

// A process to compute the mean read length of a FASTQ
process extract_libtype {
    label 'bash'
    tag "$id"
   
    input:
        tuple val(id), path(samlmon_json)

    output:
        tuple val(id), env(LIBTYPE), emit: tuple_id_libtype
   
    script:
        """
            LIBTYPE=\$(grep expected_format ${samlmon_json} | awk '{print \$2}' | tr -d '",\n')
        """

}

/*
 * Transform A to G bases in FASTQ reads for hyper-editing detection
 * Converts all adenosine (A) nucleotides to guanosine (G) in sequence lines
 * Handles both single-end and paired-end reads
 */
process transform_bases_fastq {
    label 'bash'
    tag "${meta.id}"
    publishDir("${output_dir}", mode:"copy", pattern: "*_AtoG.fastq.gz")
    
    input:
        tuple val(meta), path(fastq)  // Can be single file or multiple files (R1, R2, singles)
        val output_dir

    output:
        tuple val(meta), path("*_AtoG.fastq.gz"), emit: converted_fastq

    script:
        """
        # Process all FASTQ files (handles both single-end and paired-end)
        for fq in ${fastq}; do
            base=\$(basename "\${fq}" .fastq.gz)
            base=\$(basename "\${base}" .fq.gz)
            
            # Decompress, convert A to G in sequence lines (every 4th line starting from line 2),
            # then recompress
            zcat "\${fq}" | awk 'NR%4==2 {gsub(/A/,"G"); gsub(/a/,"g")} {print}' | gzip > "\${base}_AtoG.fastq.gz"
        done
        """
}

/*
 * Transform A to G bases in reference genome FASTA for hyper-editing detection
 * Converts all adenosine (A) nucleotides to guanosine (G) in sequence lines (non-header lines)
 */
process transform_bases_fasta {
    label 'bash'
    tag "${fasta.baseName}"
    publishDir("${output_dir}", mode:"copy", pattern: "*_AtoG.fa")
    
    input:
        path fasta
        val output_dir

    output:
        path "*_AtoG.fa", emit: converted_fasta

    script:
        def is_compressed = fasta.name.endsWith('.gz')
        def base_name = fasta.baseName.replaceAll(/\.fa(sta)?/, '')
        def output_name = "${base_name}_AtoG.fa"
        
        if (is_compressed) {
            """
            # Decompress, convert A to G in non-header lines, save uncompressed
            zcat ${fasta} | awk '/^>/ {print; next} {gsub(/A/,"G"); gsub(/a/,"g"); print}' > ${output_name}
            """
        } else {
            """
            # Convert A to G in non-header lines (lines not starting with >)
            awk '/^>/ {print; next} {gsub(/A/,"G"); gsub(/a/,"g"); print}' ${fasta} > ${output_name}
            """
        }
}