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
    tag "drip"
    publishDir("${params.outdir}/drip/${prefix}", mode:"copy", pattern: "*.tsv")
    
    input:
        val(meta_tsv)
        val prefix

    output:
        path("*_AG.tsv"), emit: editing_ag
        path("*.tsv"), emit: editing_all

    script:
        def list = meta_tsv
        def args = []
        
        // Process list by pairs: [meta, file, meta, file, ...]
        for (int i = 0; i < list.size(); i += 2) {
            def m = list[i]  // Don't reuse 'meta' variable name
            def file = list[i + 1]
            def group = m.group ?: "group_unknown"
            def sample = m.sample ?: "sample_unknown"
            def replicate = m.rep ?: "rep1"  // Note: it's 'rep' not 'replicate' in meta
            args.add("${file}:${group}:${sample}:${replicate}")
        }
        
        def args_str = args.join(" ")

        """
        drip.py drip_${prefix} ${args_str}
        """
}