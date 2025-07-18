process pluviometer {
    label "pluviometer"
    publishDir("${params.outdir}/feature_edits", mode: "copy")
    tag "${meta.id}"

    input:
        tuple(val(meta), path(site_edits))
        path(gff)
        val(tool_format)

    output:
        tuple(val(meta), path("features.tsv"), path("aggregates.tsv"), path("pluviometer.log"), emit: tuple_sample_feature_edits)


    script:
        base_name = site_edits.BaseName
        """
        pluviometer.py \
            --sites ${site_edits} \
            --gff ${gff} \
            --format ${tool_format} \
            --cov 1 \
            --edit_threshold ${params.edit_threshold} \
            --threads ${task.cpus} \
            --aggregation_mode ${params.aggregation_mode}
        """
}