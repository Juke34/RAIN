process drip {
    label "drip"
    tag "drip_${prefix}_${tool}"
    publishDir("${params.outdir}/drip/${prefix}", mode:"copy", pattern: "*/*")
    
    input:
        tuple(val(tool), val(meta_tsv))
        val prefix
        val samples_min
        val samples_pct
        val group_samples_pct
        val group_samples
        val group_samples_edited
        val group_samples_pct_edited
        val bps

    output:
        path("*_espr/*.tsv"), emit: editing_all_espr
        path("*_espf/*.tsv"), emit: editing_all_espf

    script:
        def list = meta_tsv
        def args = []
        def countAwareArgs = params.barometer_stat_test == "beta-binomial" ? "--preserve-covered-features" : ""
        
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
        def samplesMinArg = samples_min != null ? "--min-samples ${samples_min}" : ""
        def samplesPctArg = samples_pct != null ? "--min-samples-pct ${samples_pct}" : ""
        def groupSamplesArg = group_samples != null ? "--min-group-samples ${group_samples}" : ""
        def groupSamplesPctArg = group_samples_pct != null ? "--min-group-samples-pct ${group_samples_pct}" : ""
        def groupSamplesEditedArg = group_samples_edited != null ? "--min-group-samples-edited ${group_samples_edited}" : ""
        def groupSamplesPctEditedArg = group_samples_pct_edited != null ? "--min-group-samples-pct-edited ${group_samples_pct_edited}" : ""
        def bpsArg = bps ? "--bps ${bps.join(',')}" : ""

        """ 
        drip.py --threads ${task.cpus} ${samplesMinArg} ${samplesPctArg} ${groupSamplesPctArg} ${groupSamplesArg} ${groupSamplesEditedArg} ${groupSamplesPctEditedArg} ${bpsArg} ${countAwareArgs} --output drip_${prefix} ${args_str} 
        """
}