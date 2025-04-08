process convert_sp_gxf2gxf {
    label 'agat'
    tag "$sample_id"
    publishDir "${params.outdir}/bamutil_clipoverlap", mode: 'copy'

    input:
        path(gxf)

    output:
        path ("*.gff"), emit: gff

    script:
        gff_file = genome_fasta.baseName.replaceAll(/\..+(\.gz)?$/, '')
        """
        agat_convert_sp_gxf2gxf.pl --gxf ${gxf} -o ${gff_file}.gff3
        """

}
