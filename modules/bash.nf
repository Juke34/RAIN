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
    tag "${meta.uid}"
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

/*
 * Create CSV file for AliNe input from converted FASTQ reads
 * Generates a CSV with columns: sample,fastq_1,fastq_2,strandedness,read_type
 */
process create_aline_csv_he {
    label 'bash'
    tag "${file_id}"

    input:
        tuple val(meta), path(fastq) 

    output:
        path "*.csv", emit: csv

    script:
        def sample_id = meta.sample_id
        file_id = meta.file_id[0] 
        def strandedness = meta.strandedness ? meta.strandedness : "auto"
        def read_type = meta.read_type
        
        if (fastq[1]) {
            // Paired-end
            """
            fastq0=\$(readlink -f ${fastq[0]})
            fastq1=\$(readlink -f ${fastq[1]})
            echo "${sample_id},\${fastq0},\${fastq1},${strandedness},${read_type}" > ${file_id}.csv
            """
        } else {
            // Single-end
            """
            fastq0=\$(readlink -f ${fastq[0]})
            echo "${sample_id},\${fastq0},,${strandedness},${read_type}" > ${file_id}.csv
            """
        }
}

/*
 * Collect CSV files for AliNe input from converted FASTQ reads
 * Generates a single CSV with columns: sample,fastq_1,fastq_2,strandedness,read_type
 */
process collect_aline_csv_he {
    label 'bash'
    publishDir("${output_dir}", mode:"copy", pattern: "*.csv")
    
    input:
        val all_csv  // List of tuples (meta, fastq_files)
        val output_dir

    output:
        path "*.csv", emit: csv

    script:

        def list_csv = []
        list_csv = all_csv
        list_csv_bash = list_csv.join(" "); // remove bracket and replace comma by space to be processed by bash
        """
        echo "sample,fastq_1,fastq_2,strandedness,read_type" > hyper_editing_samples.csv
        for entry in ${list_csv_bash}; do
            cat \$entry >> hyper_editing_samples.csv
        done
        """
}