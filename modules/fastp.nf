process fastp {
    label 'fastp'
    tag "${meta.uid}"
    publishDir "${params.output}/${meta.uid}/qc", mode: 'copy'
    
    input:
        tuple val(meta), path(illumina)
        val(phred_type)
    
    output:
        tuple val(meta), path("*_R?_clean.fastq")
        path("${meta.uid}_fastp_report.html")
    
    script:
        """
        if [ "${phred_type}" == "64" ]
        then
          fastp -i ${illumina[0]} -I ${illumina[1]} \
              -o ${meta.uid}_R1_clean.fastq -O ${meta.uid}_R2_clean.fastq \
              --phred64 --thread ${task.cpu} \
              --detect_adapter_for_pe --html ${meta.uid}_fastp_report.html
        else
          fastp -i ${illumina[0]} -I ${illumina[1]} \
              -o ${meta.uid}_R1_clean.fastq -O ${meta.uid}_R2_clean.fastq \
              --detect_adapter_for_pe --html ${meta.uid}_fastp_report.html --thread ${task.cpu}
        fi
        """
}

process fastp_hybrid {
    label 'fastp'
    tag "${meta.uid}"
    publishDir "${params.output}/${id}/qc", mode: 'copy'

    input:
        tuple val(meta), path(illuminaR1), path(illuminaR2), path(ont)
        val(phred_type)

    output:
        tuple val(meta), path("${meta.uid}_R1_clean.fastq"),path("${meta.uid}_R2_clean.fastq"), path(ont), emit: trimmed_hybrid
        path("${meta.uid}_fastp_report.html")

    script:
        """
        if [ "${phred_type}" == "64" ]
        then
          fastp -i ${illuminaR1} -I ${illuminaR2} \
              -o ${meta.uid}_R1_clean.fastq -O ${meta.uid}_R2_clean.fastq \
              --phred64 --thread ${task.cpu} \
              --detect_adapter_for_pe --html ${meta.uid}_fastp_report.html
        else
          fastp -i ${illuminaR1} -I ${illuminaR2} \
              -o ${meta.uid}_R1_clean.fastq -O ${meta.uid}_R2_clean.fastq \
              --detect_adapter_for_pe --html ${meta.uid}_fastp_report.html --thread ${task.cpu}
        fi
        """
}
