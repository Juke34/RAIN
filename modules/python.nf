/*
Here are described all processes related to rust
*/

/*
 * Restore original read sequences in BAM files from FASTQ
 * Used after alignment with A-to-G converted sequences to restore original bases
 */
process restore_original_sequences {
    label "pluviometer"
    tag "${bam.baseName}"
    publishDir("${output_dir}", mode:"copy", pattern: "*_restored.bam")
    
    input:
        tuple( val(meta), path(bam), path(bam_unmapped))
        val output_dir

    output:
        tuple val(meta), path("*_restored.bam"), emit: restored_bam

    script:
        def output_name = "${bam.baseName}_restored.bam"
        
        """       
        restore_sequences.py -b ${bam} -u ${bam_unmapped} -o ${output_name}
        """
}

process drip {
    label "pluviometer"
    tag "drip_${tool}"
    publishDir("${params.outdir}/drip/${prefix}", mode:"copy", pattern: "*/*")
    
    input:
        tuple(val(tool), val(meta_tsv))
        val prefix
        val samples_pct
        val group_pct

    output:
        path("*_espr/*.tsv"), emit: editing_all_espr
        path("*_espf/*.tsv"), emit: editing_all_espf

    script:
        def list = meta_tsv
        def args = []
        
        // Process list of [meta, file] pairs from groupTuple
        list.each { pair ->
            def m = pair[0]  // meta dictionary
            def file = pair[1]  // file path
            def group = m.group ?: "group_unknown"
            def sample = m.sample ?: "sample_unknown"
            def replicate = m.rep ?: "rep1"
            args.add("${file}:${group}:${sample}:${replicate}")
        }
        
        def args_str = args.join(" ")

        """ 
        drip.py --threads ${task.cpus} --min-samples-pct ${samples_pct} --min-group-pct ${group_pct} --output drip_${prefix} ${args_str} 
        """
}