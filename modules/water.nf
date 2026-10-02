process pluviometer {
    label "pluviometer"
    publishDir("${params.outdir}/pluviometer/${tool_format}/", mode: "copy")
    tag "${meta.uid}_${tool_format}"

    input:
        tuple(val(meta), val(tool_format), path(site_edits))
        path(gff)

    output:
        tuple(val(meta), val(tool_format), path("*features.tsv"), emit: tuple_sample_feature)
        tuple(val(meta), val(tool_format), path("*aggregates.tsv"), emit: tuple_sample_aggregate)
        tuple(val(meta), val(tool_format), path("*pluviometer.log"), emit: tuple_sample_log)

    script:
        base_name = site_edits.BaseName
        """    
        pluviometer_wrapper.py \
            --sites ${site_edits} \
            --gff ${gff} \
            --format ${tool_format} \
            --cov ${params.cov_threshold} \
            --edit_threshold ${params.edit_threshold} \
            --threads ${task.cpus} \
            --aggregation_mode ${params.aggregation_mode} \
            --output "${meta.uid}_${tool_format}"
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

/*
 * Barometer - Exhaustive biomarker analysis of RAIN editing data.
 *
 * Consumes the drip aggregates + features TSV files (per edit type, e.g. AG, AC,
 * and per value type espf/espr) and runs barometer_analyze.py (differential /
 * multivariate / ML analysis) followed by barometer_report.py (interactive HTML
 * report). One invocation per (editType, valueType) pair.
 */
process barometer_analyze {
    label "pluviometer"
    tag "${editType}_${valueType}"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), val(valueType), path(aggregates), path(features)

    output:
        tuple val(editType), val(valueType), path("barometer_results/*"), emit: results
        path("barometer_*.log"), emit: log

    script:
        """
        barometer_analyze.py \\
            -a ${aggregates} \\
            -f ${features} \\
            -o barometer_results \\
            -j ${task.cpus} \\
            --stat-test ${params.barometer_stat_test} \\
            --max-bmks ${params.barometer_max_bmks} \\
            &> barometer_analyze_${valueType}.log
        """
}

/*
 * Barometer - Exhaustive biomarker analysis of RAIN editing data.
 *
 * Consumes the drip aggregates + features TSV files (per edit type, e.g. AG, AC,
 * and per value type espf/espr) and runs barometer_report.py (interactive HTML
 * report). One invocation per (editType).
 */
process barometer_report {
    label "pluviometer"
    tag "${editType}"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), path(espf), path(espr)

    output:
        path("barometer_report.html"), emit: report
        path("barometer_*.log"), emit: log

    script:
        """
        barometer_report.py \\
            --espf ${espf} \\
            --espr ${espr} \\
            -o barometer_report.html \\
            --embed-images \\
            &> barometer_report.log
        """
}

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