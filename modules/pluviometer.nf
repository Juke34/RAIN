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
