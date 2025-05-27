/*
Here are described all processes related to bash
*/

// A process to compute the mean read length of a FASTQ
process extract_libtype {
    label 'bash'
    tag "$id"
   
    input:
        tuple val(id), path(samlmon_json)

    output:
        tuple val(id), env(LIBTYPE), emit: tuple_id_libtype
   
    script:
        """
            LIBTYPE=\$(grep expected_format ${samlmon_json} | awk '{print \$2}' | tr -d '",\n')
        """

}