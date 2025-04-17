process normalize_gxf {
    label 'agat'
    publishDir "${params.outdir}/agat_gff3", mode: 'copy'

    input:
        path(gxf)

    output:
        path ("*.gff3"), emit: gff

    script:
        base_name = gxf.baseName.replaceAll(/\..+(\.gz)?$/, '')
        """
        agat config --expose --tabix
        agat_convert_sp_gxf2gxf.pl --gxf ${gxf} -o ${base_name}_normalized.gff3
        """
}
