process drip {
    label "drip"
    tag "drip_${prefix}_${tool}"
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

        """ 
        drip.py --threads ${task.cpus} --min-samples-pct ${samples_pct} --min-group-pct ${group_pct} ${countAwareArgs} --output drip_${prefix} ${args_str} 
        """
}