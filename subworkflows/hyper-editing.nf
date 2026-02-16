/*
 * RAIN Subworkflow: Hyper-editing detection
 * 
 * Detects extensive RNA editing events through A-to-G conversion strategy
 * This approach captures heavily edited reads that may escape standard detection
 */

include { AliNe as ALIGNMENT } from "${baseDir}/modules/aline.nf"
include { convert_to_fastq; samtools_fasta_index; samtools_split_mapped_unmapped } from "${baseDir}/modules/samtools.nf"
include { transform_bases_fastq; transform_bases_fasta } from "${baseDir}/modules/bash.nf"
include { multiqc } from "${baseDir}/modules/multiqc.nf"
include { restore_original_sequences } from "${baseDir}/modules/python.nf"

/**
 * Hyper-editing discovery workflow
 * 
 * Strategy:
 * 1. Extract unmapped reads from primary alignment
 * 2. Convert adenosines to guanosines in reads and reference
 * 3. Re-align with converted sequences
 * 4. Restore original bases for downstream analysis
 * 
 */
workflow HYPER_EDITING {
    
    take:
        unmapped_bam_chunks    // Unmapped read chunks from primary alignment
        genome        // Genomic reference sequence
        aline_profile // AliNe profile in coma-separated format
        clean_annotation // Annotation file for AliNe
        quality_threshold       // Quality score filter threshold
        output_he         // output directory path  ier
    
    main:
        // For now, create an empty channel until processes are implemented
        // This prevents the workflow from failing while we build out the functionality
        Channel.empty().set { filtered_bams }
        
        // Stage 1: Convert unmapped BAM to FASTQ format
        extracted_reads = convert_to_fastq(unmapped_bam_chunks, output_he)
        
        // Stage 2: Apply A-to-G nucleotide transformation to reads
        converted_reads = transform_bases_fastq(extracted_reads, output_he)
        
        // Stage 3: Generate A-to-G converted reference genome
        converted_reference = transform_bases_fasta(genome, output_he)
       
        // Stage 4: Collect converted reads paths for AliNe
        converted_reads
            .map { meta, fastq -> fastq }
            .collect()
            .map { files -> 
                // Create a comma-separated list or single file path
                files.collect { it.toString() }.join(',')
            }
            .set { reads_path }

        // Stage 5: Build alignment index for converted reference
        alignment_index = samtools_fasta_index(converted_reference)

        // Stage 6: Perform alignment with converted sequences
        ALIGNMENT (
            'Juke34/AliNe -r v1.6.0', // Select pipeline
            "${workflow.resume?'-resume':''} -profile ${aline_profile}", // workflow opts supplied as params for flexibility
            "-config ${params.aline_profiles}",
            reads_path,
            converted_reference,
            "--read_type ${params.read_type}",
            "--aligner ${params.aligner}",
            "--strandedness ${params.strandedness}",
            clean_annotation,
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
            .set { aligned_bams }  // Store the channel

        // create a channel of tuples with (meta, bam, fastq) joined by the id:chr21_small_R1 and sample name   
        aligned_bams.map { id, bam -> tuple(id, bam) }
                    .join(
                        extracted_reads.map { meta2, fastq -> tuple(meta2.id, meta2, fastq) }
                    )
                    .map { id, bam, meta2, fastq -> tuple(meta2, bam, fastq) }
                    .set { sample_bam_fastq } 

        // Stage 8: Reconstruct original sequences from converted alignments
        reconstructed_bams = restore_original_sequences(sample_bam_fastq, output_he)
       
        // Stage 9: Split BAM into mapped and unmapped reads
        samtools_split_mapped_unmapped(reconstructed_bams)
    
    emit:
        bam_unmapped = samtools_split_mapped_unmapped.out.unmapped_bam
        bam_mapped = samtools_split_mapped_unmapped.out.mapped_bam
}
