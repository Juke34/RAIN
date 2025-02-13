process multiqc {
    label 'multiqc'
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path log_files
    path multiqc_config

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    multiqc -p . -c ${multiqc_config}
    """
}