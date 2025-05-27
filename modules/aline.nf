/*
Adapted from
from https://github.com/mahesh-panchal/nf-cascade
*/
process AliNe {
    tag "$pipeline_name"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val pipeline_name     // String
        val profile           // String
        val config
        val reads
        val genome
        val read_type
        val aligner
        val library_type
        val cache_dir          // String

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
                reads,
                "--reference ${genome}",
                read_type,
                aligner,
                library_type,
                "--outdir $task.workDir/AliNe",
        ].join(" ")
        // Copy command to shell script in work dir for reference/debugging.
        file("$task.workDir/nf-cmd.sh").text = nxf_cmd.join(" ")
        // Run nextflow command locally
        def process = nxf_cmd.execute(null, cache_path.toFile())
        process.waitFor()
        stdout = process.text
        assert process.exitValue() == 0: stdout
        // Copy nextflow log to work directory
        cache_path.resolve(".nextflow.log").copyTo("${task.workDir}/nextflow.log")

    output:
        path "AliNe"  , emit: output
        val stdout, emit: log
}