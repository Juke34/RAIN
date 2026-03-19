/*
Adapted from
from https://github.com/mahesh-panchal/nf-cascade
*/
process AliNe {
    tag "$pipeline_name"
    label 'aline'
    publishDir "${output_dir}", mode: 'copy'
    
    errorStrategy 'terminate' // to avoid any retry
    maxRetries 0  // Override global retry config - do not retry this process

    input:
        val pipeline_name     // String
        val profile           // String
        val config
        val reads
        val genome
        val read_type
        val aligner
        val library_type
        val annotation
        val cache_dir          // cache directory
        val output_dir        // output directory 

    when:
        task.ext.when == null || task.ext.when

    exec:
        def cache_path = file(cache_dir)
        assert cache_path.mkdirs()
        // construct nextflow command
        def nxf_cmd = [
            'nextflow run',
                pipeline_name,
                profile,
                config,
                "--reads ${reads}",
                "--reference ${genome}",
                read_type,
                aligner,
                library_type,
                "--annotation ${annotation}",
                "--data_type rna",
                "--outdir $task.workDir/AliNe",
        ].join(" ")
        // Copy command to shell script in work dir for reference/debugging.
        file("$task.workDir/nf-cmd.sh").text = nxf_cmd
        // Run nextflow command locally in cache directory
        def process = nxf_cmd.execute(null, cache_path.toFile())
        // Print process output to stdout and stderr
        process.consumeProcessOutput(System.out, System.err)
        process.waitFor()
        stdout = process.text
        // Copy nextflow log to work directory
        cache_path.resolve(".nextflow.log").copyTo("${task.workDir}/nextflow.log")
        assert process.exitValue() == 0: stdout

    output:
        path "AliNe"  , emit: output
        val stdout, emit: log
}
