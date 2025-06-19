#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

// Input/output params
params.reads        = null // "/path/to/reads_{1,2}.fastq.gz/or/folder"
params.genome       = null // "/path/to/genome.fa"
params.annotation   = null // "/path/to/annotations.gff3"
params.outdir       = "rain_result"
params.clipoverlap  = false

/* Specific AliNe params (some are shared with RAIN)*/

// Read feature params
read_type_allowed        = [ 'short_paired', 'short_single', 'pacbio', 'ont' ]
params.read_type         = null // short_paired, short_single, pacbio, ont
strandedness_allowed     = [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ] // see https://github.com/Juke34/AliNe for more information
params.strandedness      = null
params.read_length       = null // Use by star to set the sjdbOverhang parameter

// Edit counting params
edit_site_tools = ["reditools2", "reditools3", "jacusa2", "sapin"]
params.edit_site_tool = "reditools3"
params.edit_threshold = 1
params.aggregation_mode = "all"

// Aline profiles
aline_profile_allowed = [ 'docker', 'singularity', 'local', 'itrop' ]

// Aline ressource config used
params.aline_profiles = "$baseDir/config/resources/custom_aline.config" // e.g. "docker, singularity,itrop,local"

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
    --reads                     path to the reads file, folder or csv. If a folder is provided, all the files with proper extension in the folder will be used. You can provide remote files (commma separated list).
                                    file extension expected : <.fastq.gz>, <.fq.gz>, <.fastq>, <.fq> or <.bam>. 
                                                              for paired reads extra <_R1_001> or <_R2_001> is expected where <R> and <_001> are optional. e.g. <sample_id_1.fastq.gz>, <sample_id_R1.fastq.gz>, <sample_id_R1_001.fastq.gz>)
                                    csv input expects 6 columns: sample, fastq_1, fastq_2, strandedness and read_type. 
                                    fastq_2 is optional and can be empty. Strandedness, read_type expects same values as corresponding AliNe parameter; If a value is provided via AliNe paramter, it will override the value in the csv file.
                                    Example of csv file:
                                        sample,fastq_1,fastq_2,strandedness,read_type
                                        control1,path/to/data1.fastq.bam,,auto,short_single
                                        control2,path/to/data2_R1.fastq.gz,path/to/data2_R2.fastq.gz,auto,short_paired
    --genome                    Path to the reference genome in FASTA format.
    --read_type                 Type of reads among this list ${read_type_allowed} (no default)

        Output:
    --output                    Path to the output directory (default: $params.outdir)

       Optional input:
    --aligner                   Aligner to use [default: $params.aligner]
    --edit_site_tool            Tool used for detecting edited sites. Default: $params.edit_site_tool
    --strandedness              Set the strandedness for all your input reads (default: null). In auto mode salmon will guess the library type for each fastq sample. [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]
    --edit_threshold            Minimal number of edited reads to count a site as edited (default: 1)
    --aggregation_mode          Mode for aggregating edition counts mapped on genomic features. See documentation for details. Options are: "all" (default) or "cds_longest"
    --clipoverlap               Clip overlapping sequences in read pairs to avoid double counting. (default: false)

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
    outdir                     : ${params.outdir}

Alignment Parameters
 aline_profiles                : ${params.aline_profiles}
    aligner                    : ${params.aligner}
    

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
include {fasta_unzip} from "$baseDir/modules/pigz.nf"
include {samtools_index; samtools_fasta_index; samtools_sort_bam} from './modules/samtools.nf'
include {reditools2} from "./modules/reditools2.nf"
include {reditools3} from "./modules/reditools3.nf"
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
      workflow.containerEngine == "singularity" ||
      workflow.containerEngine == "docker"
  ) { "executer selected" }
else { exit 1, "No executer selected: please use a profile activating docker or singularity (e.g. -profile docker/singularity/itrop)"}

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

// ----------------------------------------------------------------------------
//                               DEAL WITH CSV FILE FIRST
// ----------------------------------------------------------------------------
        def path_reads     = params.reads
        def via_url = false
        def via_csv = false
        Channel.empty().set{csv_ch}

        if ( path_reads.endsWith('.csv') ){
            log.info "Using CSV input file: ${path_reads}"
            //   --------- BAM CSV  CASE ---------
            via_csv = true
            File input_csv = new File(path_reads)
            if(!input_csv.exists()){ 
                error "The input ${path_reads} file does not exist!\n" 
            }
            
            csv_ch = Channel.fromPath(params.reads)
                                .splitCsv(header: true, sep: ',')
                                .map { row ->
                                    if(row.sample == null || row.sample.trim() == ''){ 
                                        error "The input ${params.reads} file does not contain a 'sample' column!\n" 
                                    } 
                                    def sample_id = row.sample
                                    if(row.input_1 == null || row.input_1.trim() == ''){ 
                                        error "The input ${params.reads} file does not contain a 'input_1' column!\n" 
                                    }
                                    def input_1 = file(row.input_1.trim())
                                    if( input_1.toString().endsWith('.bam') || Rain_utilities.is_fastq(input_1) ) { 
                                        
                                        def input_type = input_1.toString().endsWith('bam') ? 'bam' : 'fastq'
                                        def input_url = null 
                                        if (! Rain_utilities.is_url(input_1) ) {
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
                                        def meta = [ id: sample_id, strandedness: libtype, input_type: input_type, is_url: input_url ]
                                        return tuple(meta, input_1)
                                    }
                                    else {
                                        error "The input ${row.input_1} file is not a BAM or FASTQ file!\n"
                                    }
                                }
        // If we are here it is because the csv is correct.
        // lets create another channel with meta information
        Channel.fromPath(params.reads)
            .splitCsv(header: true)
            .map { row ->
                    def File file = new File(row.input_1)
                    fileClean = file.baseName.replaceAll(/\.(gz)$/, '') // remove .gz
                    fileClean = fileClean.replaceAll(/\.(fastq|fq)$/, '') // remove .fastq or .fq

                return tuple (fileClean, row.sample, row.strandedness)
                }
            .set { csv_meta_ch }
        }

        // Separate FASTQ samples to BAM samples
        csv_in_bam   = csv_ch.filter { meta, reads -> meta.input_type == 'bam' }
        csv_in_fastq = csv_ch.filter { meta, reads -> meta.input_type != 'fastq' }

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
                            if (  Rain_utilities.is_url(bam_path_reads) ) {
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
                bams = my_samples.flatten().map { it -> [it.name.split('_')[0], it] }
                                .groupTuple()
                                .ifEmpty { exit 1, "Cannot find reads matching ${bam_path_reads}!\n" }
            } else {
                if (bam_path_reads.endsWith('.bam')) {
                    bams = Channel.fromFilePairs(bam_path_reads, size: 1, checkIfExists: false)
                }
            }

            // Set structure with dictionary as first value
            bams = bams.map {  id, bam_file -> 
                        def strand = params.strandedness
                        strand = (strand && strand.toUpperCase() != "AUTO") ? strand : null // if strandedness is set to auto, set it to null

                        def meta = [ id: id, strandedness: strand ]
                        tuple(meta, bam_file)
                    } 
        }

        // sort the bam files
        Channel.empty().set{sorted_bam}
        sorted_bam = samtools_sort_bam( bams )

// ----------------------------------------------------------------------------
//                               DEAL WITH FASTQ FILES
// ----------------------------------------------------------------------------
        // Perform AliNe alignment
        Channel.empty().set{aline_alignments_all}
        def fastq_list=[]
        def aline_data_in = null

        if ( via_csv ){
            aline_data_in = path_reads
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
                        if (  Rain_utilities.is_fastq( it ) ){
                            if (  Rain_utilities.is_url(it) ) {
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
                        if ( Rain_utilities.is_fastq( path_reads ) ){
                            log.info "The input ${path_reads} is a fastq file!"
                            if (  Rain_utilities.is_url(path_reads) ) {
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
                'Juke34/AliNe -r v1.5.1', // Select pipeline
                "${workflow.resume?'-resume':''} -profile ${aline_profile}", // workflow opts supplied as params for flexibility
                "-config ${params.aline_profiles}",
                "--reads ${aline_data_in}",
                genome,
                "--read_type ${params.read_type}",
                "--aligner ${params.aligner}",
                "--strandedness ${params.strandedness}",
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
            
            // SET CORRECT ID NAME WHEN CSV catched from CSV sample row
            if (via_csv){
                csv_meta_ch
                    .join(aline_alignments)
                    .map{ name_AliNe, sample, strandedness, bam ->
                            return tuple (sample, bam)
                        }
                    .set { aline_alignments }
            }

            // GET TUPLE [ID, OUTPUT_SALMON_LIBTYPE] FILES
            if ( params.strandedness && params.strandedness.toUpperCase() != "AUTO" ) {
                log.info "Strandedness type is set to ${params.strandedness}, no need to extract it from salmon output"
                aline_alignments_all = aline_alignments.map { id, bam -> 
                                                                def meta = [ id: id, strandedness: params.strandedness ]
                                                                            tuple(meta, bam)}
            } else {
                log.info "Try to get strandedness from AliNe salmon output"
                aline_libtype = ALIGNMENT.out.output
                    .map { dir ->
                        files("$dir/salmon_strandedness/*/*.json", checkIfExists: true)  // Find BAM files inside the output directory
                    }
                    .flatten()  // Ensure we emit each file separately
                    .map { json -> 
                                def name = json.getParent().getName().split('_seqkit')[0]  // Extract the base name of the BAM file. _seqkit is the separator. The name is in the fodler containing the json file. Why take this one? Because it is the same as teh bam name set by Aline. It will be used to sync both values
                                tuple(name, json)
                        }  // Convert each BAM file into a tuple, with the base name as the first element

                // Extract the library type from the JSON file
                aline_libtype = extract_libtype(aline_libtype)
                
                // because id and name can be different in csv case, we add the aline name to the tuple
                aline_alignments_tmp = aline_alignments.map { id, file -> 
                            def aline_name = file.getName().split('_seqkit')[0]  // Extract the base name of the BAM file. _seqkit is the separator.
                            tuple(aline_name, id, file)
                        }

                aline_alignments_tmp.join(aline_libtype, remainder: true)
                    .map { aline_name, id, file, lib -> 
                            def meta = [ id: id, strandedness: lib ]
                            tuple(meta, file)
                        }  // Join the two channels on the key (the name of the BAM file)
                    .set { aline_alignments_all }  // Store the channel
                
                // Here strandedness is null if it was not guessed by AliNe. Try to catch it in CSV if csv case
                if (via_csv) {
                    sample_with_strand = aline_alignments_all.filter { meta, reads -> meta.strandedness && meta.strandedness != 'auto' }
                    sample_no_strand   = aline_alignments_all.filter { meta, reads -> !meta.strandedness || meta.strandedness == 'auto' }

                    log.info "Try to get strandedness from CSV file"
                    // Get the strandedness from the csv file
                    csv_meta_ch.map { read_id, sample_id, strandedness -> 
                        [sample_id, strandedness] 
                    }
                    .set { csv_meta_ch_short }  // Set the channel with the new strandedness

                    channel_data_tojoin = sample_no_strand.map { meta, bam_path -> 
                                            [meta.id, [meta, bam_path]] 
                                        }

                    channel_data_tojoin
                    .join(csv_meta_ch_short)
                    .map { id, data_pair, strandedness ->
                            def (meta, bam_path) = data_pair
                            meta.strandedness = strandedness
                            [meta, bam_path]
                        }
                    .set{sample_no_strand}

                    aline_alignments_all = sample_with_strand.concat(sample_no_strand)
                }
            }
        }

        // MERGE ALINE BAM AND INPUT BAM TOGETHER
        tuple_sample_sortedbam = aline_alignments_all.mix(sorted_bam)
        log.info "The following bam file(s) will be processed by RAIN:"
        tuple_sample_sortedbam.view()

        // STEP 1 QC with fastp ?
        Channel.empty().set{logs} 
        // stat on aligned reads
        fastqc_ali(tuple_sample_sortedbam, "ali")
        logs.concat(fastqc_ali.out).set{logs} // save log
        // remove duplicates
        gatk_markduplicates(tuple_sample_sortedbam)
        logs.concat(gatk_markduplicates.out.log).set{logs} // save log
        // stat on bam without duplicates∂
        fastqc_dup(gatk_markduplicates.out.tuple_sample_dedupbam, "dup")
        logs.concat(fastqc_dup.out).set{logs} // save log
        // Clip overlap
        if (params.clipoverlap) {
            bamutil_clipoverlap(gatk_markduplicates.out.tuple_sample_dedupbam)
            tuple_sample_bam_processed = bamutil_clipoverlap.out.tuple_sample_clipoverbam
            // stat on bam with overlap clipped
            fastqc_clip(tuple_sample_bam_processed, "clip")
            logs.concat(fastqc_clip.out).set{logs} // save log
        } else {
            tuple_sample_bam_processed = gatk_markduplicates.out.tuple_sample_dedupbam
        }
        // index bam
        samtools_index(tuple_sample_bam_processed)
        // report with multiqc
        // multiqc(logs.collect(),params.multiqc_config)

        // Select site detection tool
        switch (params.edit_site_tool) {
            case "jacusa2":
                // Create a fasta index file of the reference genome
                samtools_fasta_index(genome.collect())
                jacusa2(samtools_index.out.tuple_sample_bam_bamindex, samtools_fasta_index.out.tuple_fasta_fastaindex.collect())
                break
            case "sapin":
                sapin(tuple_sample_bam_processed, genome.collect())
                break
            case "reditools2":
                reditools2(samtools_index.out.tuple_sample_bam_bamindex, genome.collect(), params.region)
                normalize_gxf(annotation.collect())
                pluviometer(reditools2.out.tuple_sample_serial_table, normalize_gxf.out.gff.collect(), "reditools2")
                break
            case "reditools3":
                reditools3(samtools_index.out.tuple_sample_bam_bamindex, genome.collect())
                normalize_gxf(annotation.collect())
                pluviometer(reditools3.out.tuple_sample_serial_table, normalize_gxf.out.gff.collect(), "reditools3")
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
