#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

// Input/output params
params.reads = "/path/to/reads_{1,2}.fastq.gz/or/folder"
params.genome = "/path/to/genome.fa"
params.outdir = "result"
params.reads_extension = ".fastq.gz" // Extension used to detect reads in folder
params.paired_reads_pattern = "_{1,2}"

// Read feature params
read_type_allowed = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type = "short_paired" // short_paired, short_single, pacbio, ont

// Aligner params
align_tools = ['hisat2']
params.aligner = 'hisat2'
params.bowtie2_options = ''
params.hisat2_options = ''
params.star_options = ''

// Report params
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

//*************************************************
// STEP 1 - LOG INFO
//*************************************************

// Header message
log.info """
IRD
.-./`) .-------.     ______
\\ .-.')|  _ _   \\   |    _ `''.
/ `-' \\| ( ' )  |   | _ | ) _  \\
 `-'`\"`|(_ o _) /   |( ''_'  ) |
 .---. | (_,_).' __ | . (_) `. |
 |   | |  |\\ \\  |  ||(_    ._) '
 |   | |  | \\ `'   /|  (_.\\.' /
 |   | |  |  \\    / |       .'
 '---' ''-'   `'-'  '-----'`


RAIN - RNA Alterations Investigation using Nextflow
===================================================

"""

if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    ********* RAIN - RNA Alterations Investigation using Nextflow *********

        Usage example:
    nextflow run main.nf --illumina short_reads_Ecoli --genus Escherichia --species coli --species_taxid 562 -profile docker -resume
    --help                      prints the help section

        Input sequences:
    --reads                     path to the illumina read file (fastq or fastq.gz) (default: $params.reads)
    --genome                    path to the genome (default: $params.genome)

        Output:
    --output                    path to the output directory (default: $params.outdir)

       Optional input:
    --aligner                Aligner to use [default: $params.aligner]
    --read_type              type of reads among this list ${read_type_allowed} (default: short_paired)
    --reads_extension        Extension of the read files [default: $params.reads_extension]

        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README [debug, test, itrop, singularity, local, docker]
    -resume                     resume the workflow where it stopped
    """
}

// Parameter message

log.info """

General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     read_type                 : ${params.read_type}
     outdir                     : ${params.outdir}
  
Alignment Parameters
 hisat2 parameters
     hisat2_options             : ${params.hisat2_options}


Report Parameters
 MultiQC parameters
     multiqc_config             : ${params.multiqc_config}

 """


//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include { AliNe as ALIGNMENT } from "./modules/aline.nf"
include {bamutil_clipoverlap} from './modules/bamutil.nf'
include {fastp} from './modules/fastp.nf'
include {fastqc as fastqc_raw; fastqc as fastqc_ali; fastqc as fastqc_dup; fastqc as fastqc_clip} from './modules/fastqc.nf'
include {hisat2_index; hisat2} from './modules/hisat2.nf'
include {gatk_markduplicates } from './modules/gatk.nf'
include {multiqc} from './modules/multiqc.nf' 
include {samtools_sam_to_bam; samtools_sort} from './modules/samtools.nf' 

//*************************************************
// STEP 3 - Deal with parameters
//*************************************************

// Check aligner params. Can be a list (comma or space separated)
def aligner_list=[]
if( !params.aligner ){
    exit 1, "Error: <aligner> parameter is empty, please provide a aligner(s) among this list ${align_tools}.\n"
} else {
    str_list = params.aligner.tokenize(',')
    str_list.each {
        str_list2 = it.tokenize(' ')
        str_list2.each {
            if ( ! (it in align_tools) ){
                exit 1, "Error: <${it}> aligner not acepted, please provide aligner(s) among this list ${align_tools}.\n"
            }
            else{
                aligner_list.add(it)
            }
        }
    }
}

// check profile
if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker')
  ) { "executer selected" }
else { exit 1, "No executer selected: -profile docker/singularity"}


//*************************************************
// STEP 4 -  Workflow
//*************************************************

workflow {
        main:

        ALIGNMENT (
            'Juke34/AliNe -r v1.1.0',         // Select pipeline
            "-profile ${workflow.profile}",   // workflow opts supplied as params for flexibility
            "--reads ${params.reads}",
            "--genome ${params.genome}",
            "--read_type ${params.read_type}",
            "--aligner ${params.aligner}",
        )
        ALIGNMENT.out.output // The results folder
        .map { dir -> file( dir.resolve("alignment_results/alignment/hisat2/*.bam"), checkIfExists: true ) } // The relative path to the sample sheet from `results/`
        .set { aline_alignments } // Name the channel
        aline_alignments.view()
        
        rain(aline_alignments)
}

workflow rain {

    take:
        tuple_sample_sortedbam

    main:

        // STEP 1 QC with fastp ?
        Channel.empty().set{logs}
 
        // stat on aligned reads
        fastqc_ali(tuple_sample_sortedbam, "ali")
        logs.concat(fastqc_ali.out).set{logs} // save log
        // remove duplicates
        gatk_markduplicates(tuple_sample_sortedbam)
        logs.concat(gatk_markduplicates.out.log).set{logs} // save log
        // stat on bam without duplicates
        fastqc_dup(gatk_markduplicates.out.tuple_sample_dedupbam, "dup")
        logs.concat(fastqc_dup.out).set{logs} // save log
        // Clip overlap
        bamutil_clipoverlap(gatk_markduplicates.out.tuple_sample_dedupbam)
        // stat on bam with overlap clipped
        fastqc_clip(bamutil_clipoverlap.out.tuple_sample_clipoverbam, "clip")
        logs.concat(fastqc_clip.out).set{logs} // save log
        // report with multiqc
        multiqc(logs.collect(),params.multiqc_config)
}