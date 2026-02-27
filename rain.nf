#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************
// ------------------------------------
/* ---- Params specific to RAIN ---- */
// ------------------------------------
// Input/output params
params.reads        = null // "/path/to/reads_{1,2}.fastq.gz/or/folder"
params.genome       = null // "/path/to/genome.fa"
params.annotation   = null // "/path/to/annotations.gff3"
params.outdir       = "rain_result"
params.clip_overlap  = false
params.clean_duplicate  = true

// Edit counting params
edit_site_tools = ["reditools2", "reditools3", "jacusa2", "sapin"]
params.edit_site_tool = "reditools3"
params.edit_threshold = 1
params.aggregation_mode = "all"
params.skip_hyper_editing = false // Skip hyper-editing detection
// Report params
params.multiqc_config = "$baseDir/config/multiqc_config.yaml" // MultiQC config file

/* Specific tool params */
params.region = "" // e.g. chr21 - Used to limit the analysis to a specific region by REDITOOLS2

// others
params.help = null
params.monochrome_logs = false // if true, no color in logs
params.debug = false // Enable debug output
params.use_slurm_for_aline = false // Whether to submit AliNe as a separate SLURM job when using an HPC environment

// --------------------------------------------------
/* ---- Params shared between RAIN and AliNe ---- */
// --------------------------------------------------

// Read feature params
read_type_allowed        = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type         = null // short_paired, short_single, pacbio, ont
strandedness_allowed     = [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ] // see https://github.com/Juke34/AliNe for more information
params.strandedness      = null
params.fastqc            = false

// -------------------------------------
/* ---- Params Specific to AliNe ---- */
// -------------------------------------
// The rest of params are in custom_aline.config.nf file
// Aline profiles
aline_profile_allowed = [ 'docker', 'singularity', 'local', 'itrop' ]
// Aline ressource config used
params.aline_profiles = "$baseDir/nextflow_aline.config" // e.g. "docker, singularity,itrop,local"
// made in aline but params here because it is main step
params.trimming_fastp = false
// Aligner params
align_tools = [ 'bbmap', 'bowtie', 'bowtie2', 'bwaaln', 'bwamem', 'bwamem2', 'bwasw', 'dragmap', 'graphmap2', 'hisat2', 'kallisto', 'last', 'minimap2', 'novoalign', 'nucmer', 'ngmlr', 'salmon', 'star', 'subread', 'sublong' ]
params.aligner = 'hisat2'
// AliNe version
params.aline_version = 'v1.6.2'
//*************************************************
// STEP 1 - HELP
//*************************************************

log.info header()
if (params.help) { exit 0, helpMSG() }


// Params check
// Check aligner params. Can be a list (comma or space separated)
def edit_site_tool_list=[]
if( !params.edit_site_tool ){
    exit 1, "Error: <edit_site_tool> parameter is empty, please provide a aligner(s) among this list ${align_tools}.\n"
} else {
    str_list = params.edit_site_tool.tokenize(',')
    str_list.each {
        str_list2 = it.tokenize(' ')
        str_list2.each {
            if ( ! (it.toLowerCase() in edit_site_tools) ){
                exit 1, "Error: <${it}> tool not acepted, please provide editing site tool among this list ${align_tools}.\n"
            }
            else{
                edit_site_tool_list.add(it.toLowerCase())
            }
        }
    }
}

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
    --genome                    Path to the reference genome in FASTA format.
    --read_type                 Type of reads among this list ${read_type_allowed} [no default]
    --reads                     path to the reads file, folder or csv. If a folder is provided, all the files with proper extension in the folder will be used. You can provide remote files (commma separated list).
                                    file extension expected : <.fastq.gz>, <.fq.gz>, <.fastq>, <.fq> or <.bam>. 
                                                              for paired reads extra <_R1_001> or <_R2_001> is expected where <R> and <_001> are optional. e.g. <sample_1.fastq.gz>, <sample_R1.fastq.gz>, <sample_R1_001.fastq.gz>)
                                    CSV Input Format:
                                        4 required columns: group, input_1, strandedness and read_type.
                                        3 optional columns: - input_2   : in case of paired-end data
                                                            - sample    : by default, the sample name is extracted from the filename and would be uniq. However, you must provide this column if you have biological or technical replicates, so that replicates from the same sample are correctly grouped. 
                                                            - replicate : Specifies the replicate number. Otherwise, each sample will be treated as independent and assigned as rep1. Required the sample column.
                                    Strandedness, read_type expects same values as corresponding AliNe parameter; If a value is provided via AliNe paramter, it will override the value in the csv file.
                                    Example of csv file:
                                        sample,input_1,input_2,strandedness,read_type,replicate
                                        control1,path/to/data1.fastq.bam,,auto,short_single,replicate1
                                        control2,path/to/data2_R1.fastq.gz,path/to/data2_R2.fastq.gz,auto,short_paired,replicate1

        Output:
    --output                    Path to the output directory [default: $params.outdir]

       Optional input:
    --aggregation_mode          Mode for aggregating edition counts mapped on genomic features. See documentation for details. Options are: "all" (default) or "cds_longest"
    --aligner                   Aligner to use [default: $params.aligner]
    --clean_duplicate           Remove PCR duplicates from BAM files using GATK MarkDuplicates. [default: $params.clean_duplicate]
    --clip_overlap              Clip overlapping sequences in read pairs to avoid double counting. [default: $params.clipoverlap]
    --debug                     Enable debug output for troubleshooting. [default: $params.debug]
    --edit_site_tool            Tool used for detecting edited sites. [default: $params.edit_site_tool]
    --edit_threshold            Minimal number of edited reads to count a site as edited [default: $params.edit_threshold]
    --fastqc                    run fastqc on main steps [default: $params.fastqc]
    --skip_hyper_editing        Skip hyper-editing detection step for unmapped reads. [default: $params.skip_hyper_editing]
    --strandedness              Set the strandedness for all your input reads [default: $params.strandedness]. In auto mode salmon will guess the library type for each fastq sample. [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]

        Nextflow options:
    -profile                    Change the profile of nextflow both the engine and executor more details on github README [debug, test, itrop, singularity, local, docker]
    """
}

// Parameter message
log.info """

General Parameters
    annotation                 : ${params.annotation}
    reads                      : ${params.reads}
    genome                     : ${params.genome}
    strandedness               : ${params.strandedness}
    read_type                  : ${params.read_type}
    hyper-editing              : ${params.skip_hyper_editing ? "skipped" : "performed"}
    clip_overlap               : ${params.clip_overlap}
    clean_duplicate            : ${params.clean_duplicate}
    outdir                     : ${params.outdir}

Alignment Parameters
    aline_config               : ${params.aline_profiles}
    aligner                    : ${params.aligner}
    

Edited Site Detection Parameters
    edit_site_tool             : ${params.edit_site_tool}
    edit_threshold             : ${params.edit_threshold}
    region                     : ${params.region} 

Report Parameters
 MultiQC parameters
     multiqc_config            : ${params.multiqc_config}

 """


//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include { AliNe as ALIGNMENT } from "./modules/aline.nf"
include {normalize_gxf} from "./modules/agat.nf"
include { extract_libtype; recreate_csv_with_abs_paths; collect_aline_csv; filter_drip_by_aggregation_mode; filter_drip_features_by_type} from "./modules/bash.nf"
include {bamutil_clipoverlap} from './modules/bamutil.nf'
include {fastp} from './modules/fastp.nf'
include {fastqc as fastqc_ali; fastqc as fastqc_dup; fastqc as fastqc_clip} from './modules/fastqc.nf'
include {gatk_markduplicates } from './modules/gatk.nf'
include {jacusa2} from "./modules/jacusa2.nf"
include {multiqc} from './modules/multiqc.nf'
include {fasta_unzip} from "$baseDir/modules/pigz.nf"
include {samtools_index; samtools_fasta_index; samtools_sort_bam as samtools_sort_bam_raw; samtools_sort_bam as samtools_sort_bam_merged; samtools_split_mapped_unmapped; samtools_merge_bams} from './modules/samtools.nf'
include {reditools2} from "./modules/reditools2.nf"
include {reditools3} from "./modules/reditools3.nf"
include {pluviometer as pluviometer_jacusa2; pluviometer as pluviometer_reditools2; pluviometer as pluviometer_reditools3; pluviometer as pluviometer_sapin} from "./modules/pluviometer.nf"
include {drip_aggregates; drip_features} from "./modules/python.nf"
include {sapin} from "./modules/sapin.nf"

include {HYPER_EDITING} from "./subworkflows/hyper-editing.nf"


//*************************************************
// STEP 3 - Deal with parameters
//*************************************************

// check RAIN profile - /!\ profile must be sync with AliNe profile as much as possible
if (
      workflow.containerEngine == "singularity" ||
      workflow.containerEngine == "docker"
  ) { "executer selected" }
else { exit 1, "No executer selected: please use a profile activating docker or singularity (e.g. -profile docker/singularity/itrop)"}

// check AliNE profile
def aline_profile_list=[]
def use_slurm_for_aline = params.use_slurm_for_aline
str_list = workflow.profile.tokenize(',')
str_list.each {
    if ( it in aline_profile_allowed ){
         aline_profile_list.add(it)
    }
}
def aline_profile = aline_profile_list.join(',')

// Check csv file
def via_csv     = false
def path_reads  = params.reads
if ( path_reads.endsWith('.csv') ){  
    via_csv = true
    println  "Using CSV input file: ${path_reads}"
}

// check read type parameter
def read_typep = params.read_type == "null" ? null : params.read_type // if read_type is "null" it will be set to null, otherwise it will be the value of params.read_type

if( via_csv ) {
    if( read_typep ){
        if ( ! (read_typep.toLowerCase() in read_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${read_typep}> read_type not acepted, please provide a read type among this list ${read_type_allowed}."
        }
        println """    This value will replace any read_type value found in your csv!"""
    } else {
        println """    No read_type provided by --read_type parameter, value will be taken from the csv file."""
    }
} else {
    if( ! read_typep ){
        exit 1, "Error: <read_type> parameter is empty it is allowed only when CSV provided as input. Please provide a read type among this list ${read_type_allowed}."
    } else {
        if ( ! (read_typep.toLowerCase() in read_type_allowed*.toLowerCase()) ){
            exit 1, "Error: <${read_typep}> read_type not acepted, please provide a read type among this list ${read_type_allowed}."
        }
    }
}

//*************************************************
// STEP 4 -  Workflow
//*************************************************

workflow {
        main:
        
        Channel.empty().set{logs} // logs channel
// ----------------------------------------------------------------------------
        // --- DEAL WITH REFERENCE ---
        // check if reference exists
        if (!params.reads) {
            exit 1, "No input reads provided with --reads (fastq or bam files as file, file list, folder or csv)."
        }
        if (!params.genome) {
            exit 1, "You must provide a path to the genome with --genome"
        }
        Channel.fromPath(params.genome, checkIfExists: true)
               .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
               .set{genome_raw}
        // unzip it if needed
        fasta_unzip(genome_raw)
        fasta_unzip.out.genomeFa.set{genome} // set genome to the output of fasta_unzip
// ----------------------------------------------------------------------------
        // --- DEAL WITH ANNOTATION ---
        Channel.empty().set{annotation}
        if (params.annotation){
        Channel.fromPath(params.annotation, checkIfExists: true)
               .ifEmpty { exit 1, "Cannot find annotation matching ${params.annotation}!\n" }
               .set{annotation}
        }
        // normalize the annotation
        normalize_gxf(annotation)
        normalize_gxf.out.gff.set{clean_annotation}
// ----------------------------------------------------------------------------
//                               DEAL WITH CSV FILE FIRST
// ----------------------------------------------------------------------------
        def via_url = false
        Channel.empty().set{csv_ch}

        if ( via_csv ){
            //   --------- BAM / FASTQ CSV CASE ---------
            File input_csv = new File(path_reads)
            if(!input_csv.exists()){ 
                error "The input ${path_reads} file does not exist!\n" 
            }
            
            csv_ch = Channel.fromPath(params.reads)
                                .splitCsv(header: true, sep: ',')
                                .map { row ->
                                    if(row.group == null || row.group.trim() == ''){ 
                                        error "The input ${params.reads} file does not contain a 'group' column!\n" 
                                    }
                                    def group = row.group.trim()
                                    if(row.input_1 == null || row.input_1.trim() == ''){ 
                                        error "The input ${params.reads} file does not contain a 'input_1' column!\n" 
                                    }
                                    def input_1 =file(RainUtils.getAbsolutePath(row.input_1))
                                    if( input_1.toString().endsWith('.bam') || RainUtils.is_fastq(input_1) ) { 
                                        
                                        def input_type = input_1.toString().endsWith('bam') ? 'bam' : 'fastq'
                                        def input_url = null 
                                        if (! RainUtils.is_url(input_1) ) {
                                            if (! input_1.exists() ) {
                                                error "The input ${input_1} file does not does not exits!\n"
                                            }
                                        } else {
                                            log.info "This bam input is an URL: ${input_1}"
                                            input_url = true
                                        }

                                        // check if strandedness is provided. Priority to params.strandedness. If not provided neither in csv then send an error
                                        def libtype
                                        if ( params.strandedness ) {
                                            libtype = params.strandedness
                                        } else {
                                            if(row.strandedness == null || row.strandedness.trim() == ''){ 
                                                error "The input ${params.reads} file does not contain a 'strandedness' column!\n" 
                                            }
                                        }
                                        libtype = row.strandedness.trim()

                                        // read type
                                        def read_type = null
                                        def pair = false
                                        if ( !read_typep ) {
                                            if (row.read_type != null) {
                                                read_type_value = row.read_type.trim().toLowerCase()
                                                if (read_type_value){
                                                    if ( ! ( read_type_value in read_type_allowed*.toLowerCase()) ){
                                                        error "The input ${input_csv} file contains an invalid read type value: ${read_type_value}. Please provide one of the following values: ${read_type_allowed}."
                                                    } else {
                                                        read_type = read_type_value
                                                    }
                                                } else {
                                                    error "The input ${input_csv} file contains an empty read_type value for input_1 ${input_1}!"
                                                }
                                            } else {
                                                error """Error: The input file ${input_csv} does not contain a read_type column, and the --read_type parameter was not provided.
    Please specify the read type either by including a read_type column in the input file or by using the --read_type option."""
                                            }
                                        } else {
                                            read_type = read_typep
                                        }
                                        // check its is paired or not
                                        def input_2
                                        if(row.input_2) {
                                            input_2 = file(RainUtils.getAbsolutePath(row.input_2))
                                        }
                                        if ( input_2 ) {
                                            if (read_type == "short_paired") {
                                                pair = true
                                                file_id2 = RainUtils.cleanPrefix(input_2)
                                            } else {
                                                log.info "The input ${input_csv} file contains a second fastq file for input_1 ${input_1} but the read_type is set to <${read_type}>! R2 will not be taken into account! paired set to false."
                                            }
                                        } else {
                                            if (read_type == "short_paired" && input_type == "fastq") {
                                                error "The input ${input_csv} file does not contain a second fastq file for input_1 ${input_1} but the read_type is set to <short_paired>!"
                                            }
                                        }
                                        // Define file_id based on the filename and read type using the utility function
                                        def file_id1 = RainUtils.cleanPrefix(input_1)
                                        // uid is similar to file_id[0] but file_id needed to make difference between R1 and R2 in paired end case.
                                        def uid = pair ? RainUtils.get_file_uid(input_1) : RainUtils.cleanPrefix(input_1)
                                        def file_id = pair ? [file_id1, file_id2] : [file_id1]
                                        def file_abspath = pair ? [input_1, input_2] : [input_1]
                                        
                                        // Deal with sample and replicate columns
                                        def sample
                                        def replicate
                                        if(row.sample == null || row.sample.trim() == ''){ 
                                            if(row.replicate == null || row.replicate.trim() == ''){
                                                // case no sample column and no replicate column, we consider each file as a sample and assign the sample as the uid
                                                log.info "No sample column and no replicate column, we consider each file as a sample and assign the sample as <${uid}>!"
                                                sample = uid
                                                replicate = "rep1"
                                            } else {
                                                error "You provided a replicate column but the sample column is missing. Please add a sample column!"
                                            } 
                                        } else {
                                            sample = row.sample.trim()
                                            if(row.replicate == null && row.replicate.trim() == ''){
                                                replicate = "rep1"
                                            } else {
                                                replicate = row.replicate.trim()
                                            }
                                        }

                                        // create meta dictionary with all the information about the sample, it will be used in the rest of the pipeline to keep trace of the different parameters for each sample
                                        def meta = [ uid: uid, file_id: file_id, sample: sample, rep: replicate, group: group, file_abspath: file_abspath, strandedness: libtype, input_type: input_type, is_url: input_url, read_type: read_type]
                                        
                                        
                                        return tuple(meta, input_1)
                                    }
                                    else {
                                        error "The input ${row.input_1} file is not a BAM or FASTQ file!\n"
                                    }
                                }
        }

        // Separate FASTQ samples to BAM samples
        csv_in_bam   = csv_ch.filter { meta, reads -> meta.input_type == 'bam' }
        csv_in_fastq = csv_ch.filter { meta, reads -> meta.input_type != 'fastq' }

        // Check read_type is provided when not using CSV
        if (!via_csv && !params.read_type) {
            exit 1, "Error: --read_type parameter is required when not using a CSV input file. Please provide a read type among ${read_type_allowed}.\n"
        }

// ----------------------------------------------------------------------------
//                               DEAL WITH BAM FILES
// ----------------------------------------------------------------------------

        def bam_path_reads = params.reads
        def bam_list=[]
        Channel.empty().set{bams}
        
        if (via_csv) {
            bams = csv_in_bam
        }
        else {

            //   --------- BAM LIST CASE ---------
            if( bam_path_reads.indexOf(',') >= 0) {
                // Cut into list with coma separator
                str_list = bam_path_reads.tokenize(',')
                // loop over elements
                str_list.each {
                    str_list2 = it.tokenize(' ')
                    str_list2.each {
                        if ( it.endsWith('.bam') ){
                            if (  RainUtils.is_url(it) ) {
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
                def File input_reads = new File(bam_path_reads)
                if(input_reads.exists()){
                    if ( input_reads.isDirectory()) {
                        log.info "The input ${bam_path_reads} is a folder! Check for bam presence..."
                        // in case of folder provided, add a trailing slash if missing
                        def matchingFiles = input_reads.listFiles().findAll { 
                            it.name ==~ /.*\.bam$/
                        }
                        if (matchingFiles) {
                            bam_path_reads = "${input_reads}" + "/*.bam"
                            log.info "Bam file found"
                        } else {
                            bam_path_reads = null
                            log.info "No bam file found"
                        }
                    }
                    //   --------- BAM FILE  CASE ---------
                    else {
                        if ( bam_path_reads.endsWith('.bam') ){
                            log.info "The input ${bam_path_reads} is a bam file!"
                            if (  RainUtils.is_url(bam_path_reads) ) {
                                log.info "This bam input is an URL: ${bam_path_reads}"
                                via_url = true
                            }
                        }
                    }
                }
            }
        }

        // Load Bam
        if (! via_csv) {
            if ( params.strandedness && params.strandedness.toUpperCase() == "AUTO" ) {
                log.info "⚠️ The <auto> Strandedness mode cannot be invoked for BAM files, it will be set to null, otherwise you need to use the --strandedness parameter. See help for more information."
            }
            if (via_url ){
                my_samples = Channel.of(bam_list)
                bams = my_samples.flatten().map { it -> 
                                                    uid = RainUtils.get_file_uid(it)
                                                    [uid, it] 
                                                }
                                .groupTuple()
                                .ifEmpty { exit 1, "Cannot find reads matching ${bam_path_reads}!\n" }
            } else {
                if (bam_path_reads.endsWith('.bam')) {
                    bams = Channel.fromFilePairs(bam_path_reads, size: 1, checkIfExists: false)
                }
            }

            // Set structure with dictionary as first value
            bams = bams.map {  uid, bam_file -> 
                        def strand = params.strandedness
                        strand = (strand && strand.toUpperCase() != "AUTO") ? strand : null // if strandedness is set to auto, set it to null
                        def meta = [ uid: uid, strandedness: strand, read_type: params.read_type ]
                        tuple(meta, bam_file)
                    } 
        }

        // stat on aligned reads
        if(params.fastqc){
            fastqc_ali(bams, "ali")
            logs.concat(fastqc_ali.out).set{logs} // save log
        }

// ----------------------------------------------------------------------------
//                               DEAL WITH FASTQ FILES
// ----------------------------------------------------------------------------
        // Perform AliNe alignment
        Channel.empty().set{aline_alignments_all}
        def fastq_list=[]
        def aline_data_in = null

        if ( via_csv ){
            aline_data_in_ch = recreate_csv_with_abs_paths(csv_ch)
            aline_data_in_ch = collect_aline_csv(aline_data_in_ch.collect(), "AliNe")
            aline_data_in = aline_data_in_ch
        } 
        else {
            //   --------- FASTQ LIST CASE ---------
            if( path_reads.indexOf(',') >= 0) {
                // Cut into list with coma separator
                str_list = path_reads.tokenize(',')
                // loop over elements
                str_list.each {
                    str_list2 = it.tokenize(' ')
                    str_list2.each {
                        if (  RainUtils.is_fastq( it ) ){
                            if (  RainUtils.is_url(it) ) {
                                log.info "This fastq input is an URL: ${it}"
                                via_url = true
                            }
                            fastq_list.add(file(it)) // use file insted of File for URL
                        }
                    }
                }
                aline_data_in = fastq_list.join(',') // Join the list into a string for AliNe input (filtered to contain only fastq files)
            }
            else {
                //   --------- FASTQ FOLDER CASE ---------
                File input_reads = new File(path_reads)
                if(input_reads.exists()){
                    if ( input_reads.isDirectory()) {
                        log.info "The input ${path_reads} is a folder! Check for fastq presence..."
                        // in case of folder provided, add a trailing slash if missing
                        def matchingFiles = input_reads.listFiles().findAll { 
                            it.name ==~ /.*_(R?1|R?2)(_.*)?\.(fastq|fq)(\.gz)?$/
                        }
                        if (matchingFiles) {
                            aline_data_in =   "${input_reads}" + "/" // Aline will handle the folder
                            log.info "Fastq file found"
                        } else {
                            log.info "No fastq file found"
                        }
                    }
                    //   --------- FASTQ FILE  CASE ---------
                    else {
                        if ( RainUtils.is_fastq( path_reads ) ){
                            log.info "The input ${path_reads} is a fastq file!"
                            if (  RainUtils.is_url(path_reads) ) {
                                log.info "This fastq input is an URL: ${path_reads}"
                                via_url = true
                            }
                            aline_data_in = path_reads // Use the file as is
                        }
                    }
                }
            }
        }
        
        Channel.empty().set{aline_alignments_all}
        if (aline_data_in){

            ALIGNMENT (
                "Juke34/AliNe -r ${params.aline_version}", // Select pipeline
                "${workflow.resume?'-resume':''} -profile ${aline_profile}", // workflow opts supplied as params for flexibility
                "-config ${params.aline_profiles}",
                aline_data_in,
                genome,
                "--read_type ${params.read_type}",
                "--aligner ${params.aligner}",
                "--strandedness ${params.strandedness}",
                clean_annotation,
                workflow.workDir.resolve('Juke34/AliNe').toUriString(),
                use_slurm_for_aline  // Pass info about whether to use slurm submission
            )

            // GET TUPLE [ID, BAM] FILES
            ALIGNMENT.out.output
                .map { dir ->
                    files("$dir/alignment/**/*.bam", checkIfExists: true)  // Find BAM files recursively inside the output directory
                }
                .flatten()  // Ensure we emit each file separately
                .map { bam -> 
                            def name = bam.getName().split('_AliNe')[0]  // Extract the base name of the BAM file. _seqkit is the separator.
                            def meta = [ uid: name ]
                            if (! via_csv){
                                meta.strandedness = params.strandedness
                                meta.read_type = params.read_type
                            }
                            tuple(meta, bam)
                    }  // Convert each BAM file into a tuple, with the base name as the first element
                .set { aline_alignments }  // Store the channel

            if (via_csv) {
                aline_alignments.map { meta, bam -> tuple(meta.uid, bam) }  // Reorder the tuple to have uid as the first element for joining
                    .join(  
                        csv_ch.map { meta, file -> tuple(meta.uid, meta) }  // Create a channel of tuples with (uid, meta) from the CSV metadata
                    )
                    .map { uid, bam, meta -> tuple(meta, bam) }  
                    .set { aline_alignments }  // Store the channel
            }

            // UPDATE STRANDEDNESS INFORMATION
            log.info "Try to get strandedness from AliNe salmon output when strandedness is set to auto or not provided"
            aline_libtype = ALIGNMENT.out.output
                .map { dir ->
                    files("$dir/salmon_strandedness/*.json")  // Find BAM files inside the output directory
                }
                .flatten()  // Ensure we emit each file separately
                .map { json -> 
                            def name = json.getName().split('_AliNe')[0]  // Extract the base name of the BAM file. _seqkit is the separator. The name is in the fodler containing the json file. Why take this one? Because it is the same as teh bam name set by Aline. It will be used to sync both values
                            tuple(name, json)
                    }  // Convert each BAM file into a tuple, with the base name as the first element

            // Extract the library type from the JSON file
            aline_libtype = extract_libtype(aline_libtype)

            // add the libtype to the meta information
            aline_alignments.map { meta, bam -> tuple(meta.uid, meta, bam) }
                .join(aline_libtype, remainder: true)
                .map { uid, meta, bam, lib ->
                        if (lib != null) {
                            meta.strandedness = lib
                        } 
                        tuple(meta, bam)
                    }  // Join the two channels on the key (the name of the BAM file)
                .set { aline_alignments }  // Store the channel
        }

        // MERGE ALINE BAM AND INPUT BAM TOGETHER
        all_input_bam = aline_alignments.mix(bams)
        if (params.debug) {
            all_input_bam.view { meta, bam -> log.info "RAIN processing ${bam} (sample: ${meta.uid})" }
        }

// -------------------------------------------------------
// ----------------- DETECT HYPER_EDITING ----------------
// -------------------------------------------------------

        // Split BAM into mapped and unmapped reads and Add _AliNe to bam file name to look similar as fastq mapped via AliNe
        samtools_split_mapped_unmapped(all_input_bam)

        // Process hyper-editing if not skipped
        Channel.empty().set{all_bam}
        if (!params.skip_hyper_editing) {
            
            HYPER_EDITING(
                samtools_split_mapped_unmapped.out.unmapped_bam,
                genome,
                aline_profile,
                use_slurm_for_aline,
                clean_annotation,
                30,  // quality threshold
                "${params.outdir}/hyper_editing",
            )
            
            // Access the results
            hyperedit_bam_mapped = HYPER_EDITING.out.bam_mapped
            bam_unmapped = HYPER_EDITING.out.bam_unmapped
         
            // create a channel of tuples with (meta, bam, bam_he) joined by the id
            samtools_split_mapped_unmapped.out.mapped_bam.map { meta, bam -> tuple(meta.uid, meta, bam) }
                        .join(
                            hyperedit_bam_mapped.map { meta2, bam_he -> tuple(meta2.uid, meta2, bam_he) }
                        )
                        .map { id, meta, bam, meta2, bam_he -> tuple(meta, bam, bam_he) }
                        .set { meta_bam_bamhe } 
            
            // Merge the final bam with the original bam in case of hyper-editing to keep all reads for edition site detection
            all_bam = samtools_merge_bams(meta_bam_bamhe, "bam_appended_with_he")
            
            // Add hyper-editing sample to analysis
            // Here we have bam containing normal and hyper_editing reads and bam containing only hyper-editing reads.
            // The bam with only hyper-editing reads will be used to count hyper-editing sites specifically. 
            //all_bam = merged_bams.mix(hyperedit_bam_mapped)
       
        } else {
            all_bam = samtools_split_mapped_unmapped.out.mapped_bam
        }
        
        // Sort the bam files
        all_bam_sorted = samtools_sort_bam_merged(all_bam)

        // remove duplicates
        if(params.clean_duplicate){
            gatk_markduplicates(all_bam_sorted)
            all_bam_sorted_dedup = gatk_markduplicates.out.tuple_sample_dedupbam
            logs.concat(gatk_markduplicates.out.log).set{logs} // save log
            // stat on bam without duplicates
            if(params.fastqc){
                fastqc_dup(all_bam_sorted_dedup, "dup")
                logs.concat(fastqc_dup.out).set{logs} // save log
            }
        }
        else {
            all_bam_sorted_dedup = all_bam_sorted
        }

        // Clip overlap
        if (params.clip_overlap) {
            bamutil_clipoverlap(gatk_markduplicates.out.tuple_sample_dedupbam)
            all_bam_sorted_dedup_clip = bamutil_clipoverlap.out.tuple_sample_clipoverbam
            // stat on bam with overlap clipped
            if(params.fastqc){
                fastqc_clip(all_bam_sorted_dedup_clip, "clip")
                logs.concat(fastqc_clip.out).set{logs} // save log
            }
        } else {
            all_bam_sorted_dedup_clip = gatk_markduplicates.out.tuple_sample_dedupbam
        }
        
        // index mapped bam
        samtools_index(all_bam_sorted_dedup_clip)
        final_bam_for_editing = samtools_index.out.tuple_sample_bam_bamindex

// -------------------------------------------------------
// ----------------- DETECT EDITING SITES ----------------
// -------------------------------------------------------

        // Select site detection tool
        if ( "jacusa2" in edit_site_tool_list ){ 
                // Create a fasta index file of the reference genome
                samtools_fasta_index(genome.collect())
                jacusa2(final_bam_for_editing, samtools_fasta_index.out.tuple_fasta_fastaindex.collect())
                pluviometer_jacusa2(jacusa2.out.tuple_sample_jacusa2_table, clean_annotation.collect(), "jacusa2")
        }
        if ( "sapin" in edit_site_tool_list ){ 
                sapin(tuple_sample_bam_processed, genome.collect())
        }
        if ( "reditools2" in edit_site_tool_list ){ 
                reditools2(final_bam_for_editing, genome.collect(), params.region)
                pluviometer_reditools2(reditools2.out.tuple_sample_serial_table, clean_annotation.collect(), "reditools2")
        }
        if ( "reditools3" in edit_site_tool_list ){ 
                reditools3(final_bam_for_editing, genome.collect())
                pluviometer_reditools3(reditools3.out.tuple_sample_serial_table, clean_annotation.collect(), "reditools3")
                if(via_csv){
                    // Deal with AGGREGATES
                    drip_aggregates(pluviometer_reditools3.out.tuple_sample_aggregate.collect())
                    filter_drip_by_aggregation_mode(drip_aggregates.out.editing_ag, "AG")
                    // Deal with FEATURES
                    drip_features(pluviometer_reditools3.out.tuple_sample_feature.collect())
                    filter_drip_features_by_type(drip_features.out.editing_ag, "AG")
                }
        }

        // ------------------- MULTIQC -----------------
        multiqc(logs.collect(),params.multiqc_config)
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
    Exit Status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    =======================================================
    ${c_reset}
    """
}
