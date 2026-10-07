/*
 * Barometer - Exhaustive biomarker analysis of RAIN editing data.
 *
 * Consumes the drip aggregates + features TSV files (per edit type, e.g. AG, AC,
 * and per value type espf/espr) and runs barometer_analyze.py (differential /
 * multivariate / ML analysis) followed by barometer_report.py (interactive HTML
 * report). One invocation per (editType, valueType) pair.
 */
process barometer_analyze {
    label "barometer"
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
    label "barometer"
    tag "${editType}"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), path(espf), path(espr), val(site)

    output:
        path("barometer_report.html"), emit: report
        path("barometer_*.log"), emit: log

    script:
        """
        barometer_report.py \\
            --espf ${espf} \\
            --espr ${espr} \\
            ${site ? "--site ${site}" : ""} \\
            -o barometer_report.html \\
            --embed-images \\
            &> barometer_report.log
        """
}

/*
 * Barometer on per-site ESPR matrices (beta-binomial); one invocation per edit type.
 * Input: standard drip.py ESPR TSVs computed on the per-site pluviometer output.
 */
process barometer_analyze_sites {
    label "barometer"
    tag "${editType}_sites"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), path(sites)

    output:
        tuple val(editType), path("barometer_results/*"), emit: results
        path("barometer_*.log"), emit: log

    script:
        """
        barometer_analyze.py \\
            --sites ${sites} \\
            -o barometer_results \\
            -j ${task.cpus} \\
            --stat-test ${params.barometer_stat_test} \\
            --max-bmks ${params.barometer_max_bmks} \\
            &> barometer_analyze_sites.log
        """
}