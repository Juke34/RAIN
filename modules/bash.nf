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
    tag "${meta.uid}"

    input:
        tuple val(meta), path(fastq) 

    output:
        path "*.csv", emit: csv

    script:
        def sample_id = meta.sample_id
        def strandedness = meta.strandedness ? meta.strandedness : "auto"
        def read_type = meta.read_type
        
        if (fastq[1]) {
            // Paired-end
            """
            fastq0=\$(readlink -f ${fastq[0]})
            fastq1=\$(readlink -f ${fastq[1]})
            echo "${sample_id},\${fastq0},\${fastq1},${strandedness},${read_type}" > ${meta.uid}.csv
            """
        } else {
            // Single-end
            """
            fastq0=\$(readlink -f ${fastq[0]})
            echo "${sample_id},\${fastq0},,${strandedness},${read_type}" > ${meta.uid}.csv
            """
        }
}

/*
 * Collect CSV files for AliNe input from converted FASTQ reads
 * Generates a single CSV with columns: sample,fastq_1,fastq_2,strandedness,read_type
 */
process collect_aline_csv {
    label 'bash'
    publishDir("${params.outdir}/${output_dir}", mode:"copy", pattern: "*.csv")
    
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
        echo "sample,fastq_1,fastq_2,strandedness,read_type" > aline_input.csv
        for entry in ${list_csv_bash}; do
            cat \$entry >> aline_input.csv
        done
        """
}

/*
 * Recreate CSV file with absolute paths from metadata
 * Takes input samples with metadata and generates a new CSV with absolute file paths
 * This allows conversion of relative paths to absolute paths for CI/CD compatibility
 */
process recreate_csv_with_abs_paths {
    label 'bash'
    tag "${meta.uid}"
    
    input:
        tuple val(meta), path(file) 

    output:
        path "*.csv", emit: csv

    script:
        
        def sample_id = meta.sample_id
        def strandedness = meta.strandedness
        def read_type = meta.read_type
        def file1 = meta.file_abspath[0]
        def file2 = meta.file_abspath[1] ? meta.file_abspath[1] : ""
        """
        echo "${sample_id},${file1},${file2},${strandedness},${read_type}" > ${meta.uid}.csv
        """
}

/*
 * Filter drip output files by AggregationMode
 * Creates separate files for each unique AggregationMode value
 * Output files are named with the AggregationMode as suffix (e.g., drip_AG_all_isoforms.tsv)
 */
process filter_drip_by_aggregation_mode {
    label 'bash'
    tag "${tsv.baseName}"
    publishDir("${params.outdir}/drip/aggregates_${output_dir}_by_mode", mode:"copy", pattern: "*.tsv")
    
    input:
        path tsv
        val output_dir

    output:
        path "*.tsv", emit: filtered_files
        val output_dir

    script:
        def basename = tsv.baseName
        """
        # Separate the 3 kinds of values with header written only once per file
        awk 'BEGIN {FS=OFS="\\t"; h1=0; h3=0}
             NR==1 {header=\$0; next}
             {
                 if(\$1=="." && \$2=="."){
                     if(h1==0){print header > "total.tsv"; h1=1}
                     print \$0 > "total.tsv"
                 } 
                 else if(\$2=="."){
                     fname="bySeq_"\$1".tsv"
                     if(!(fname in h2)){print header > fname; h2[fname]=1}
                     print \$0 > fname
                 } 
                 else {
                     if(h3==0){print header > "therest"; h3=1}
                     print \$0 > "therest"
                 }
             }' ${tsv}
        
        # Get unique AggregationMode values (column 6)
        # Skip header, extract column 6, get unique values, remove empty lines
        MODES=\$(tail -n +2 therest | cut -f 6 | sort -u | grep -v '^\$')
        
        # Extract header
        HEADER=\$(head -n 1 ${tsv})

        # Create a file for each AggregationMode
        for MODE in \$MODES; do
            OUTPUT_FILE="${basename}_\${MODE}.tsv"
            
            # Write header
            echo "\$HEADER" > "\$OUTPUT_FILE"
            
            # Filter rows by this AggregationMode (column 6)
            tail -n +2 therest | awk -F'\t' -v mode="\$MODE" '\$6 == mode' >> "\$OUTPUT_FILE"
            
            echo "Created \$OUTPUT_FILE with \$(tail -n +2 \"\$OUTPUT_FILE\" | wc -l) rows"
        done
        
        echo "Filtering complete. Created files for modes: \$MODES"
        """
}

/*
 * Filter drip_features output files by Type
 * Creates separate files for each unique Type value
 * Output files are named with the Type as suffix (e.g., drip_features_AG_gene.tsv)
 */
process filter_drip_features_by_type {
    label 'bash'
    tag "${tsv.baseName}"
    publishDir("${params.outdir}/drip/features_${output_dir}_by_type", mode:"copy", pattern: "*_*.tsv")
    
    input:
        path tsv
        val output_dir

    output:
        path "*_*.tsv", emit: filtered_files

    script:
        def basename = tsv.baseName
        """
        # Extract header
        HEADER=\$(head -n 1 ${tsv})
        
        # Get unique Type values (column 4)
        # Skip header, extract column 4, get unique values, remove empty lines
        TYPES=\$(tail -n +2 ${tsv} | cut -f 4 | sort -u | grep -v '^\$')
        
        # Create a file for each Type
        for TYPE in \$TYPES; do
            OUTPUT_FILE="${basename}_\${TYPE}.tsv"
            
            # Write header
            echo "\$HEADER" > "\$OUTPUT_FILE"
            
            # Filter rows by this Type (column 4)
            tail -n +2 ${tsv} | awk -F'\t' -v type="\$TYPE" '\$4 == type' >> "\$OUTPUT_FILE"
            
            echo "Created \$OUTPUT_FILE with \$(tail -n +2 \"\$OUTPUT_FILE\" | wc -l) rows"
        done
        
        echo "Filtering complete. Created files for types: \$TYPES"
        """
}