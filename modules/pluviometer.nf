process pluviometer {
    label "pluviometer"
    publishDir("${params.outdir}/feature_edits", mode: "copy")
    tag "${meta.id}"

    input:
        tuple(val(meta), path(site_edits))
        path(gff)
        val(tool_format)

    output:
        tuple(val(meta), path("*.feature_edits.tsv"), emit: tuple_sample_feature_edits)


    script:
        base_name = site_edits.BaseName
        """
        # cp ${workflow.projectDir}/bin/stats/*.py ./
        python ${workflow.projectDir}/bin/stats/pluviometer.py \
            --sites ${site_edits} \
            --gff ${gff} \
            --format ${tool_format} \
            --cov 1 \
            --edit_threshold ${params.edit_threshold} \
            --aggregation_mode ${params.aggregation_mode} \
            --output ${base_name}.feature_edits.tsv
        """
}