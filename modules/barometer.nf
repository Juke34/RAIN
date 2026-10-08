/*
 * Barometer - Exhaustive biomarker analysis of RAIN editing data.
 *
 * Each DRIP output (aggregates / features / sites × espf / espr) is analysed
 * independently by a single barometer_analyze invocation. Six runs per edit
 * type, fully parallel. Results are then pooled by barometer_merge into a
 * truly-global ranking, and barometer_report renders the interactive HTML.
 */
process barometer_analyze {
    label "barometer"
    tag "${editType}_${vtype}_${mtype}"
    publishDir("${params.outdir}/barometer/${editType}/barometer_analyze/${vtype}", mode: "copy")

    input:
        tuple val(editType), val(vtype), val(mtype), path(input_file)

    output:
        tuple val(editType), val(vtype), val(mtype), path("barometer_${vtype}_${mtype}/"), emit: results
        path("barometer_*.log"), emit: log

    script:
        def flag = mtype == "aggregate" ? "-a" : (mtype == "feature" ? "-f" : "--sites")
        """
        barometer_wrapper.py \\
            ${flag} ${input_file} \\
            -o barometer_${vtype}_${mtype} \\
            -j ${task.cpus} \\
            --stat-test ${params.barometer_stat_test} \\
            --max-bmks ${params.barometer_max_bmks} \\
            &> barometer_${vtype}_${mtype}.log
        """
}

/*
 * Barometer - "truly global" cross-mtype / cross-vtype ranking.
 *
 * Consumes the per-mtype global_ranking CSVs produced by the six independent
 * barometer_analyze runs, merges them into a single global ranking, and runs
 * a second-pass analysis by redistributing the significant biomarkers back
 * to their value type. One invocation per edit type.
 *
 * The raw DRIP TSVs are passed so the second pass can rebuild the full rows.
 */
process barometer_merge {
    label "barometer"
    tag "${editType}_merge"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), path(results), path(raw_files), val(raw_keys)

    output:
        tuple val(editType), path("barometer_merged"), emit: results
        path("barometer_*.log"), emit: log

    script:
        // Local work-dir paths of the raw DRIP TSVs (Nextflow interpolates the
        // list to space-separated local paths), aligned with raw_keys.
        def localPaths = "${raw_files}".split(/\s+/)
        def rawMap = [:]
        for (i in 0..<raw_keys.size()) {
            def vt = raw_keys[i].split('/')[0]
            def mtype = raw_keys[i].split('/')[1]
            def key = mtype == "aggregate" ? "aggregates" : (mtype == "feature" ? "features" : "sites")
            rawMap[vt] = rawMap[vt] ?: [:]
            rawMap[vt][key] = localPaths[i]
        }
        // Valid JSON object: { "espf": {"aggregates": "...", ...}, "espr": {...} }
        def rawInputsJson = groovy.json.JsonOutput.toJson(rawMap)
        def resultsStr = results.collect { it.toString() }.join(" ")
        """
        barometer_wrapper.py \\
            --merge \\
            --results-dir ${resultsStr} \\
            --raw-inputs '${rawInputsJson}' \\
            -o barometer_merged \\
            --stat-test ${params.barometer_stat_test} \\
            --max-bmks ${params.barometer_max_bmks} \\
            &> barometer_merge.log
        """
}

/*
 * Barometer - Interactive HTML report.
 *
 * Consumes the per-vtype result directories (aggregates + features + sites)
 * and the merged (truly-global) ranking directory. One invocation per edit
 * type.
 */
process barometer_report {
    label "barometer"
    tag "${editType}"
    publishDir("${params.outdir}/barometer/${editType}", mode: "copy")

    input:
        tuple val(editType), path(espf), path(espr), path(merged)

    output:
        path("barometer_report.html"), emit: report
        path("barometer_*.log"), emit: log

    script:
        def espfStr = espf.collect { it.toString() }.join(" ")
        def esprStr = espr.collect { it.toString() }.join(" ")
        """
        barometer_wrapper.py \\
            --report \\
            --espf ${espfStr} \\
            --espr ${esprStr} \\
            --merged ${merged} \\
            -o barometer_report.html \\
            --embed-images \\
            &> barometer_report.log
        """
}
