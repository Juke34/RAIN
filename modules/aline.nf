/*
from https://github.com/mahesh-panchal/nf-cascade
*/
process AliNe {
    tag "$pipeline_name"

    input:
    val pipeline_name     // String
    val profile     // String
    val reads       
    val genome       
    val read_type
    val aligner

    when:
    task.ext.when == null || task.ext.when

    exec:
    // def args = task.ext.args ?: ''
    def cache_dir = java.nio.file.Paths.get(workflow.workDir.resolve(pipeline_name).toUri())
    java.nio.file.Files.createDirectories(cache_dir)
    // construct nextflow command
    def nxf_cmd = [
        'nextflow run',
            pipeline_name,
            profile,
            reads,
            genome,
            read_type,
            aligner,
            "--outdir $task.workDir/results",
    ]
    // Copy command to shell script in work dir for reference/debugging.
    file("$task.workDir/nf-cmd.sh").text = nxf_cmd.join(" ")
    // Run nextflow command locally
    def builder = new ProcessBuilder(nxf_cmd.join(" ").tokenize(" "))
    builder.directory(cache_dir.toFile())
    process = builder.start()
    assert process.waitFor() == 0: process.text
    // Copy nextflow log to work directory
    file("${cache_dir.toString()}/.nextflow.log").copyTo("$task.workDir/.nextflow.log")

    output:
    path "results"  , emit: output
    val process.text, emit: log
}