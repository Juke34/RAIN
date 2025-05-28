#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

// Input/output params
params.reads      = null // "/path/to/reads_{1,2}.fastq.gz/or/folder"
params.bam        = null // "/path/to/reads.bam/or/folder"
params.csv        = null // "/path/to/reads.bam/or/folder"
params.genome     = null // "/path/to/genome.fa"
params.annotation = null // "/path/to/annotations.gff3"
params.outdir     = "result"

/* Specific AliNe params (some are shared with RAIN)*/

// Read feature params
read_type_allowed = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type = "short_paired" // short_paired, short_single, pacbio, ont
params.library_type = "auto" // can be 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' - see https://github.com/Juke34/AliNe for more information
params.bam_library_type = null // can be 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' - see https://github.com/Juke34/AliNe for more information

// Edit counting params
edit_site_tools = ["reditools2", "jacusa2", "sapin"]
params.edit_site_tool = "reditools2"
params.edit_threshold = 1
params.aggregation_mode = "all"

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
    nextflow run rain.nf -profile docker --genome /path/to/genome.fa --annotation /path/to/annotation.gff3 --reads /path/to/reads_folder --output /path/to/output --aligner hisat2

        Parameters:
    --help                      Prints the help section

        Input sequences:
    --annotation                Path to the annotation file (GFF or GTF)
    --bam                       Path to the bam file or folder.
    --csv                       Path to the csv file containing fastq/bam files and their metadata. The csv file must contain a 'sample' column with the sample name, a 'input_1' column with the path to the bam file/or fastq1 file, a 'input_2' column with the path to R2 fastq file in case of paired data (either empty), and a 'strandedness' column with the library type (e.g. U, IU, MU, OU, ISF, ISR, MSF, MSR, OSF, OSR). 
    --reads                     Path to the reads file or folder. If a folder is provided, all the files with proper extension (.fq, .fastq, .fq.gz, .fastq.gz) in the folder will be used. You can provide remote files (commma separated list).
    --genome                    Path to the reference genome in FASTA format.

        Output:
    --output                    Path to the output directory (default: $params.outdir)

       Optional input:
    --aligner                   Aligner to use [default: $params.aligner]
    --edit_site_tool            Tool used for detecting edited sites. Default: $params.edit_site_tool
    --read_type                 Type of reads among this list ${read_type_allowed} (default: short_paired)
    --library_type              Set the library_type of your fastq reads (default: auto). In auto mode salmon will guess the library type for each sample. [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]
    --bam_library_type          Set the library_type of your BAM reads (default: null). When using BAM you must set a library type: [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]
    --edit_threshold            Minimal number of edited reads to count a site as edited (default: 1)
    --aggregation_mode          Mode for aggregating edition counts mapped on genomic features. See documentation for details. Options are: "all" (default) or "cds_longest"

        Nextflow options:
    -profile                    Change the profile of nextflow both the engine and executor more details on github README [debug, test, itrop, singularity, local, docker]
    """
}

// Parameter message

log.info """

General Parameters
    annotation                 : ${params.annotation}
    bam                        : ${params.bam}
    csv                        : ${params.csv}
    reads                      : ${params.reads}
    genome                     : ${params.genome}
     
    reads_library_type         : ${params.library_type}
    bam_library_type           : ${params.bam_library_type}
    outdir                     : ${params.outdir}

Alignment Parameters
 aline_profiles                : ${params.aline_profiles}
    aligner                    : ${params.aligner}
    reads_library_type         : ${params.library_type}

Edited Site Detection Parameters
    edit_site_tool             : ${params.edit_site_tool}
    edit_threshold             : ${params.edit_threshold}

Report Parameters
 MultiQC parameters
     multiqc_config            : ${params.multiqc_config}

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
include {samtools_index; samtools_fasta_index; samtools_sort_bam} from './modules/samtools.nf'
include {reditools2} from "./modules/reditools2.nf"
include {jacusa2} from "./modules/jacusa2.nf"
include {sapin} from "./modules/sapin.nf"
include {normalize_gxf} from "./modules/agat.nf"
include {pluviometer} from "./modules/pluviometer.nf"

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
                exit 1, "Error: <${it}> aligner not accepted, please provide aligner(s) among this list ${align_tools}.\n"
            }
            else{
                aligner_list.add(it)
            }
        }
    }
}

// Check edit site tool params. Does not accept list yet, but validates input.
if ( ! (params.edit_site_tool in edit_site_tools) ){
                exit 1, "Error: <${it}> edit site tool not accepted, please provide a tool in this list ${edit_site_tools}.\n"
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
// ----------------------------------------------------------------------------
        // --- DEAL WITH REFERENCE ---
        // check if reference exists
        Channel.fromPath(params.genome, checkIfExists: true)
               .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
               .set{genome_raw}
        // uncompress it if needed
        fasta_uncompress(genome_raw)
        fasta_uncompress.out.genomeFa.set{genome_ch} // set genome to the output of fasta_uncompress
// ----------------------------------------------------------------------------
        // --- DEAL WITH ANNOTATION ---
        Channel.empty().set{annotation_ch}
        if (params.annotation){
        Channel.fromPath(params.annotation, checkIfExists: true)
               .ifEmpty { exit 1, "Cannot find annotation matching ${params.annotation}!\n" }
               .set{annotation_ch}
        }
// ----------------------------------------------------------------------------
        def path_csv = params.csv
        def path_bam = params.bam
        def path_reads = params.reads
        def has_fastq = false
        // ---------------------- DEAL WITH BAM FILES ----------------------
        Channel.empty().set{sorted_bam}     
        if ( path_bam || path_csv ) {

            def bam_list=[]          
            def via_url = false
            def via_csv = false

            if ( path_csv && path_csv.endsWith('.csv') ){
                log.info "Using CSV input file: ${path_csv}"
                //   --------- BAM CSV  CASE ---------
                via_csv = true
                File input_csv = new File(path_csv)
                if(!input_csv.exists()){ 
                    error "The input ${params.csv} file does not exist!\n" 
                }

                bams = Channel.fromPath(params.csv)
                                    .splitCsv(header: true, sep: ',')
                                    .map { row ->
                                        if(row.sample == null || row.sample.trim() == ''){ 
                                            error "The input ${params.csv} file does not contain a 'sample' column!\n" 
                                        } 
                                        def sample_id = row.sample
                                        if(row.input_1 == null || row.input_1.trim() == ''){ 
                                            error "The input ${params.csv} file does not contain a 'input_1' column!\n" 
                                        }
                                        if( row.input_1.toString().endsWith('.bam') ){ 
                                            def input_1 = file(row.input_1.trim())
                                        
                                            if (! Rain_utilities.is_url(input_1) ) {
                                                if (! input_1.exists() ) {
                                                    error "The input ${input_1} file does not does not exits!\n"
                                                }
                                            } else {
                                                log.info "This bam input is an URL: ${input_1}"
                                            }

                                            if(row.strandedness == null || row.strandedness.trim() == ''){ 
                                                error "The input ${params.csv} file does not contain a 'strandedness' column!\n" 
                                            }
                                            def libtype      = row.strandedness ?: 'U'
                                            def meta = [ id: sample_id, libtype: libtype ]
                                            return tuple(meta, input_1)
                                        }
                                        else if ( Rain_utilities.is_fastq(row.input_1.trim()) ) {
                                            has_fastq = true
                                        /*    def input_1 = file(row.input_1.trim())
                                            if (! Rain_utilities.is_url(input_1) ) {
                                                if (! input_1.exists() ) {
                                                    error "The input ${input_1} file does not does not exits!\n"
                                                }
                                            } else {
                                                log.info "This fastq input is an URL: ${input_1}"
                                            }
                                        */
                                        }
                                        else {
                                            error "The input ${row.input_1} file is not a BAM or FASTQ file!\n"
                                        }
                                    }
            }
            else {

                //   --------- BAM LIST CASE ---------
                if( path_bam.indexOf(',') >= 0) {
                    // Cut into list with coma separator
                    str_list = path_bam.tokenize(',')
                    // loop over elements
                    str_list.each {
                        str_list2 = it.tokenize(' ')
                        str_list2.each {
                            if ( it.endsWith('.bam') ){
                                if (  Rain_utilities.is_url(it) ) {
                                    log.info "This bam input is an URL: ${it}"
                                    via_url = true
                                }
                                bam_list.add(file(it)) // use file insted of File for URL
                            }
                        }
                    }
                }
                else {
                    //   --------- BAM FOLDER CASE ---------
                    File input_reads = new File(path_bam)
                    if(input_reads.exists()){
                        if ( input_reads.isDirectory()) {
                            log.info "The input ${path_bam} is a folder!"
                            // in case of folder provided, add a trailing slash if missing
                            path_bam = "${input_reads}" + "/*.bam"
                        }
                        //   --------- BAM FILE  CASE ---------
                        else {
                            if ( path_bam.endsWith('.bam') ){
                                log.info "The input ${path_bam} is a bam file!"
                                if (  Rain_utilities.is_url(path_bam) ) {
                                    log.info "This bam input is an URL: ${path_bam}"
                                    via_url = true
                                }
                            }
                        }
                    }
                }
            }

            // Add lib
            if (! via_csv) {
                if (via_url ){
                    my_samples = Channel.of(bam_list)
                    bams = my_samples.flatten().map { it -> [it.name.split('_')[0], it] }
                                    .groupTuple()
                                    .ifEmpty { exit 1, "Cannot find reads matching ${path_reads}!\n" }
                } else {
                    bams = Channel.fromFilePairs(path_bam, size: 1, checkIfExists: true)
                        .ifEmpty { exit 1, "Cannot find reads matching ${path_bam}!\n" }
                }

                // check if bam library type is provided when bam files are used
                bams.collect().map { list -> 
                                        if (!list.isEmpty() && !params.bam_library_type ) {
                                            exit 1, "⚠️  When using BAM files, the library type must be provided. Please use the parameter --bam_library_type <libtype>.\n"
                                        }
                                    }
                // Set structure with dictionary as first value
                bams = bams.map {  id, bam_file -> 
                            def meta = [ id: id, libtype: params.bam_library_type ]
                            tuple(meta, bam_file)
                        } 
            }
            // sort the bam files
            sorted_bam = samtools_sort_bam( bams )
        }
// ----------------------------------------------------------------------------
        // DEAL WITH FASTQ FILES
        // Perform AliNe alignment
        Channel.empty().set{aline_alignments_all}
        if (path_reads || ( path_csv && has_fastq) ) {
            
            if ( path_csv && path_csv.endsWith('.csv') ){  
                path_reads = path_csv
            }

            ALIGNMENT (
                'Juke34/AliNe -r v1.4.0', // Select pipeline
                "${workflow.resume?'-resume':''} -profile ${aline_profile}", // workflow opts supplied as params for flexibility
                "-config ${params.aline_profiles}",
                "--reads ${path_reads}",
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
        all_bams = aline_alignments_all.mix(sorted_bam)
        log.info "The following bam file(s) will be processed by RAIN:"
        all_bams.view()
        rain(all_bams, genome_ch, annotation_ch)
}

workflow rain {

    take:
        tuple_sample_sortedbam
        genome
        annnotation

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
        // Create a fasta index file of the reference genome
        samtools_fasta_index(genome)

        // Select site detection tool
        switch (params.edit_site_tool) {
            case "jacusa2":
                jacusa2(samtools_index.out.tuple_sample_bam_bamindex, samtools_fasta_index.out.tuple_fasta_fastaindex)
                break
            case "sapin":
                sapin(bamutil_clipoverlap.out.tuple_sample_clipoverbam, genome)
                break
            case "reditools2":
                reditools2(samtools_index.out.tuple_sample_bam_bamindex, genome, params.region)
                normalize_gxf(annnotation)
                pluviometer(reditools2.out.tuple_sample_serial_table, normalize_gxf.out.gff, "reditools")
                break
            default:
                exit(1, "Wrong edit site tool was passed")
        }
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