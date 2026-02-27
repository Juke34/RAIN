process pluviometer {
    label "pluviometer"
    publishDir("${params.outdir}/pluviometer/${tool_format}/raw", mode: "copy")
    tag "${meta.uid}"

    input:
        tuple(val(meta), path(site_edits))
        path(gff)
        val(tool_format)

    output:
        tuple(val(meta), path("*features.tsv"), emit: tuple_sample_feature)
        tuple(val(meta), path("*aggregates.tsv"), emit: tuple_sample_aggregate)
        tuple(val(meta), path("*pluviometer.log"), emit: tuple_sample_log)

    script:
        base_name = site_edits.BaseName
        """
        pluviometer_wrapper.py \
            --sites ${site_edits} \
            --gff ${gff} \
            --format ${tool_format} \
            --cov 1 \
            --edit_threshold ${params.edit_threshold} \
            --threads ${task.cpus} \
            --aggregation_mode ${params.aggregation_mode} \
            --output "${meta.uid}_${tool_format}" 
        """
}