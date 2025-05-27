process fastp {
    label 'fastp'
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/qc", mode: 'copy'
    
    input:
        tuple val(meta), path(illumina)
        val(phred_type)
    
    output:
        tuple val(meta), path("*_R?_clean.fastq")
        path("${meta.id}_fastp_report.html")
    
    script:
        """
        if [ "${phred_type}" == "64" ]
        then
          fastp -i ${illumina[0]} -I ${illumina[1]} \
              -o ${meta.id}_R1_clean.fastq -O ${meta.id}_R2_clean.fastq \
              --phred64 \
              --detect_adapter_for_pe --html ${meta.id}_fastp_report.html
        else
          fastp -i ${illumina[0]} -I ${illumina[1]} \
              -o ${meta.id}_R1_clean.fastq -O ${meta.id}_R2_clean.fastq \
              --detect_adapter_for_pe --html ${meta.id}_fastp_report.html
        fi
        """
}

process fastp_hybrid {
    label 'fastp'
    tag "${meta.id}"
    publishDir "${params.output}/${id}/qc", mode: 'copy'

    input:
        tuple val(meta), path(illuminaR1), path(illuminaR2), path(ont)
        val(phred_type)

    output:
        tuple val(meta), path("${meta.id}_R1_clean.fastq"),path("${meta.id}_R2_clean.fastq"), path(ont), emit: trimmed_hybrid
        path("${meta.id}_fastp_report.html")

    script:
        """
        if [ "${phred_type}" == "64" ]
        then
          fastp -i ${illuminaR1} -I ${illuminaR2} \
              -o ${meta.id}_R1_clean.fastq -O ${meta.id}_R2_clean.fastq \
              --phred64 \
              --detect_adapter_for_pe --html ${meta.id}_fastp_report.html
        else
          fastp -i ${illuminaR1} -I ${illuminaR2} \
              -o ${meta.id}_R1_clean.fastq -O ${meta.id}_R2_clean.fastq \
              --detect_adapter_for_pe --html ${meta.id}_fastp_report.html
        fi
        """
}
