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


/* Specific AliNe params (some are shared with RAIN)*/

// Read feature params
read_type_allowed = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type = "short_paired" // short_paired, short_single, pacbio, ont
params.library_type = "auto" // can be 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' - see https://github.com/Juke34/AliNe for more information
// Aline profiles
aline_profile_allowed = [ 'docker', 'singularity', 'local', 'itrop' ]

// Aline ressource config used
params.aline_profiles = "$baseDir/config/ressources/custom_aline.config" // e.g. "docker, singularity,itrop,local"

// Aligner params
align_tools = ['hisat2', "STAR"]
params.aligner = 'hisat2'
params.bowtie2_options = ''
params.hisat2_options = ''
params.star_options = ''

/* Specific tool params */
params.region = "" // e.g. chr21 - Used to limit the analysis to a specific region by REDITOOLS2 

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
     annotation                 : ${params.annotation}
     reads                      : ${params.reads}
     read_type                  : ${params.read_type}
     outdir                     : ${params.outdir}
  
Alignment Parameters
 aline_profiles                 : ${params.aline_profiles}
     aligner                    : ${params.aligner}
     library_type               : ${params.library_type}
     paired_reads_pattern       : ${params.paired_reads_pattern}
     reads_extension            : ${params.reads_extension}

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
include {gatk_markduplicates } from './modules/gatk.nf'
include {multiqc} from './modules/multiqc.nf' 
include {samtools_index; samtools_fasta_index} from './modules/samtools.nf'
include {reditools2} from "./modules/reditools2.nf"
include {jacusa2} from "./modules/jacusa2.nf"
include {sapin} from "./modules/sapin.nf"

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

// check RAIN profile - /!\ profile must be sync with AliNe profile as much as possible
if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker')
  ) { "executer selected" }
else { exit 1, "No executer selected: -profile docker/singularity"}

// check AliNE profile
def aline_profile_list=[]
str_list = workflow.profile.tokenize(',')
str_list.each {
    if ( it in aline_profile_allowed ){
         aline_profile_list.add(it)
    }
}
def aline_profile = aline_profile_list.join(',')

//*************************************************
// STEP 4 -  Workflow
//*************************************************

workflow {
        main:

        ALIGNMENT (
            'Juke34/AliNe -r v1.3.0',         // Select pipeline
            "-profile ${aline_profile}",   // workflow opts supplied as params for flexibility
            "-config ${params.aline_profiles}",
            "--reads ${params.reads}",
            "--genome ${params.genome}",
            "--read_type ${params.read_type}",
            "--aligner ${params.aligner}",
            "--library_type ${params.library_type}"
        )
        ALIGNMENT.out.output
            .map { dir -> 
                files("$dir/alignment/*/*.bam", checkIfExists: true)  // Find BAM files inside the output directory
            }
            .flatten()  // Ensure we emit each file separately
            .map { bam -> tuple(bam.baseName, bam) }  // Convert each BAM file into a tuple, with the base name as the first element
            .set { aline_alignments }  // Store the channel
        
        // aline_alignments.view()
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
        // index bam
        samtools_index(bamutil_clipoverlap.out.tuple_sample_clipoverbam)
        // report with multiqc
        // multiqc(logs.collect(),params.multiqc_config)
        // Detect RNA editing with reditools2
        reditools2(samtools_index.out.tuple_sample_bam_bamindex, params.genome, params.region)
        // Create a fasta index file of the reference genome
        samtools_fasta_index(params.genome)
        jacusa2(samtools_index.out.tuple_sample_bam_bamindex, samtools_fasta_index.out.tuple_fasta_fastaindex)
        sapin(bamutil_clipoverlap.out.tuple_sample_clipoverbam, params.genome)
}