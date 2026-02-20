process normalize_gxf {
    label 'agat'
    tag "${base_name}"
    publishDir "${params.outdir}/agat_gff3", mode: 'copy'

    input:
        path(gxf)

    output:
        path ("*.gff3"), emit: gff

    script:
        base_name = RainUtils.cleanPrefix(gxf)
        """
        agat config --expose --tabix
        agat_convert_sp_gxf2gxf.pl --gxf ${gxf} -o ${base_name}_normalized.gff3
        """
}
