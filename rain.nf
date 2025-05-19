#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

// Input/output params
params.reads = null // "/path/to/reads_{1,2}.fastq.gz/or/folder"
params.bams = null // "/path/to/reads.bam/or/folder"
params.genome = "/path/to/genome.fa"
params.annotation = "/path/to/annotations.gff3"
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

// other
params.help = null
params.monochrome_logs = false // if true, no color in logs

//*************************************************
// STEP 1 - HELP
//*************************************************

log.info header()
if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    RAIN - RNA Alterations Investigation using Nextflow - v${workflow.manifest.version}

        Usage example:
    nextflow run main.nf --illumina short_reads_Ecoli --genus Escherichia --species coli --species_taxid 562 -profile docker -resume
    --help                      prints the help section

        Input sequences:
    --reads                     path to the illumina read file (fastq or fastq.gz) (default: $params.reads)
    --genome                    path to the genome (default: $params.genome)

        Annotation input:
    --annotation                path to a GFF3 file with annotations of genomic features

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
include { extract_libtype } from "./modules/bash.nf"
include {bamutil_clipoverlap} from './modules/bamutil.nf'
include {fastp} from './modules/fastp.nf'
include {fastqc as fastqc_raw; fastqc as fastqc_ali; fastqc as fastqc_dup; fastqc as fastqc_clip} from './modules/fastqc.nf'
include {gatk_markduplicates } from './modules/gatk.nf'
include {multiqc} from './modules/multiqc.nf'
include {fasta_uncompress} from "$baseDir/modules/pigz.nf"
include {samtools_index; samtools_fasta_index} from './modules/samtools.nf'
include {reditools2} from "./modules/reditools2.nf"
include {jacusa2} from "./modules/jacusa2.nf"
include {sapin} from "./modules/sapin.nf"
include {normalize_gxf} from "./modules/agat.nf"

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

        // --- DEAL WITH REFERENCE ---
        // check if reference exists
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set{genome_raw}
        // uncompress it if needed
        fasta_uncompress(genome_raw)
        fasta_uncompress.out.genomeFa.set{genome_ch} // set genome to the output of fasta_uncompress

        // DEAL WITH BAM FILES
        if ( params.bams ) {

            def bam_list=[]
            def path_bams = params.bams
            def via_URL = false
            if( path_bams.indexOf(',') >= 0) {
                // Cut into list with coma separator
                str_list = path_bams.tokenize(',')
                // loop over elements
                str_list.each {
                    str_list2 = it.tokenize(' ')
                    str_list2.each {
                        if ( it.endsWith('.bam') ){
                            if (it.startsWith('https:') || it.startsWith('s3:') || it.startsWith('az:') || it.startsWith('gs:') || it.startsWith('ftp:') ) {
                                log.info "This bam input is an URL: ${it}"
                                via_URL = true
                            }
                            bam_list.add(file(it)) // use file insted of File for URL
                        }
                    }
                }
            }
            else {
                File input_reads = new File(path_bams)
                if(input_reads.exists()){
                    if ( input_reads.isDirectory()) {
                        log.info "The input ${path_bams} is a folder!\n"
                        // in case of folder provided, add a trailing slash if missing
                        path_bams = "${input_reads}" + "/"
                    }
                    else {
                        log.info "The input ${path_bams} is a file!\n"
                        if ( path_bams.endsWith('.bam') ){
                            if (path_bams.startsWith('https:') || path_bams.startsWith('s3:') || path_bams.startsWith('az:') || path_bams.startsWith('gs:') || path_bams.startsWith('ftp:') ) {
                                log.info "This bam input is an URL: ${it}"
                                via_URL = true
                            }
                            bam_list = path_bams
                        }
                    }
                }
            }
            if (via_URL ){
                my_samples = Channel.of(bam_list)
                bams = my_samples.flatten().map { it -> [it.name.split('_')[0], it] }
                                .groupTuple()
                                .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
            } else {
                bams = Channel.fromFilePairs("${path_bams}/*.bam", size: 1, checkIfExists: true)
                    .ifEmpty { exit 1, "Cannot find reads matching ${path_bams}!\n" }
            }
            bams
        }

        // DEAL WITH FASTQ FILES
        // Perform AliNe alignment
        if (params.reads) {

            ALIGNMENT (
                'Juke34/AliNe -r v1.4.0', // Select pipeline
                "${workflow.resume?'-resume':''} -profile ${aline_profile}", // workflow opts supplied as params for flexibility
                "-config ${params.aline_profiles}",
                "--reads ${params.reads}",
                genome_ch,
                "--read_type ${params.read_type}",
                "--aligner ${params.aligner}",
                "--library_type ${params.library_type}",
                workflow.workDir.resolve('Juke34/AliNe').toUriString()
            )

            // GET TUPLE [ID, BAM] FILES
            ALIGNMENT.out.output
                .map { dir ->
                    files("$dir/alignment/*/*.bam", checkIfExists: true)  // Find BAM files inside the output directory
                }
                .flatten()  // Ensure we emit each file separately
                .map { bam -> 
                            def name = bam.getName().split('_seqkit')[0]  // Extract the base name of the BAM file. _seqkit is the separator.
                            tuple(name, bam)
                    }  // Convert each BAM file into a tuple, with the base name as the first element
                .set { aline_alignments }  // Store the channel
            
            if (params.library_type.contains("auto") ) {
                log.info "Library type is set to auto, extracting it from salmon output"
                // GET TUPLE [ID, OUTPUT_SALMON_LIBTYPE] FILES
                ALIGNMENT.out.output
                    .map { dir ->
                        files("$dir/salmon_libtype/*/*.json", checkIfExists: true)  // Find BAM files inside the output directory
                    }
                    .flatten()  // Ensure we emit each file separately
                    .map { json -> 
                                def name = json.getParent().getName().split('_seqkit')[0]  // Extract the base name of the BAM file. _seqkit is the separator. The name is in the fodler containing the json file. Why take this one? Because it is the same as teh bam name set by Aline. It will be used to sync both values
                                tuple(name, json)
                        }  // Convert each BAM file into a tuple, with the base name as the first element
                    .set { aline_libtype }  // Store the channel
                // Extract the library type from the JSON file
                aline_libtype = extract_libtype(aline_libtype)
                aline_alignments.join(aline_libtype)
                    .map { key, val1, val2 -> tuple(key, val1, val2) }
                    .set { aline_alignments_all }
            } else {
                log.info "Library type is set to ${params.library_type}, no need to extract it from salmon output"
                aline_alignments_all = aline_alignments.map { name, bam -> tuple(name, bam, params.library_type) }
            }

            // transform [ID, BAM, LIBTYPE] into [[id: 'ID', libtype: 'LIBTYPE'], file('BAM')]
            aline_alignments_all = aline_alignments_all.map { id, file, lib ->
                def meta = [ id: id, libtype: lib ]
                tuple(meta, file)
            }
        }
        // call rain
        rain(aline_alignments_all, genome_ch)
}

workflow rain {

    take:
        tuple_sample_sortedbam
        genome

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
        reditools2(samtools_index.out.tuple_sample_bam_bamindex, genome, params.region)
        // Create a fasta index file of the reference genome
        samtools_fasta_index(genome)
        jacusa2(samtools_index.out.tuple_sample_bam_bamindex, samtools_fasta_index.out.tuple_fasta_fastaindex)
        sapin(bamutil_clipoverlap.out.tuple_sample_clipoverbam, genome)
        normalize_gxf(params.annotation)
}


//*************************************************
def header(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    return """
    -${c_dim}--------------------------------------------------${c_reset}-
    ${c_blue}.-./`) ${c_white}.-------.    ${c_red} ______${c_reset}
    ${c_blue}\\ .-.')${c_white}|  _ _   \\  ${c_red} |    _ `''.${c_reset}     French National   
    ${c_blue}/ `-' \\${c_white}| ( ' )  |  ${c_red} | _ | ) _  \\${c_reset}    
    ${c_blue} `-'`\"`${c_white}|(_ o _) /  ${c_red} |( ''_'  ) |${c_reset}    Research Institute for    
    ${c_blue} .---. ${c_white}| (_,_).' __ ${c_red}| . (_) `. |${c_reset}
    ${c_blue} |   | ${c_white}|  |\\ \\  |  |${c_red}|(_    ._) '${c_reset}    Sustainable Development
    ${c_blue} |   | ${c_white}|  | \\ `'   /${c_red}|  (_.\\.' /${c_reset}
    ${c_blue} |   | ${c_white}|  |  \\    / ${c_red}|       .'${c_reset}
    ${c_blue} '---' ${c_white}''-'   `'-'  ${c_red}'-----'`${c_reset}
    ${c_purple} RAIN - RNA Alterations Investigation using Nextflow - v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

/**************         onComplete         ***************/

workflow.onComplete {

    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.success) {
        log.info "\n${c_green}    RAIN pipeline complete!"
    } else {
        log.error "${c_red}Oops .. something went wrong}"
    }

    log.info "    The results are available in the ‘${params.outdir}’ directory."
    
    def dateComplete = workflow.complete.format("dd-MMM-yyyy HH:mm:ss")
    def duration     = workflow.duration
    def succeeded    = workflow.stats['succeededCount'] ?: 0
    def cached       = workflow.stats['cachedCount'] ?: 0
    log.info """
    ================ RAIN Pipeline Summary ================
    Completed at: ${dateComplete}
    UUID        : ${workflow.sessionId}
    Duration    : ${duration}
    Succeeded   : ${succeeded}
    Cached      : ${cached}
    =======================================================
    ${c_reset}
    """
}