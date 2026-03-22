process reditools3 {
    label "reditools3"
    publishDir("${params.outdir}/reditools3", mode: "copy")
    tag "${meta.uid}"

    input:
        tuple(val(meta), path(bam), path(bamindex))
        path genome

    output:
        tuple(val(meta), path("${base_name}.site_edits_reditools3.txt"), emit: tuple_sample_serial_table)
        path("${base_name}.reditools3.log", emit: log)

    script:
        // Set the strand orientation parameter from the library type parameter
        // Terms explained in https://salmon.readthedocs.io/en/latest/library_type.html
        if (meta.strandedness in ["ISR", "SR"]) {
            // First-strand oriented
            strand_orientation = "2"
        } else if (meta.libtype in ["ISF", "SF"]) {
            // Second-strand oriented
            strand_orientation = "1"
        } else if (meta.strandedness in ["IU", "U"]) {
            // Unstranded
            strand_orientation = "0"
        } else {
            // Unsupported: Pass the library type string so that it's reported in
            // the reditools error message
            print("invalid strand ${meta.strandedness}")
            strand_orientation =0
        }
        base_name = bam.BaseName

        """
        # Trap and catch up to avoid reditools not stopping properly
        # [ERROR] (<class 'ValueError'>) MD tag not present
        # [ERROR] Killing job
        
        trap "kill 0" EXIT

        LOG=${base_name}.reditools3.log
        OUT=${base_name}.site_edits_reditools3.txt

        python -m reditools analyze ${bam} \
            --reference ${genome} \
            --strand ${strand_orientation} \
            --output-file \$OUT \
            --threads ${task.cpus} \
            --verbose \
            &> \$LOG  || exit 1

        # 🔥 Détection erreurs silencieuses REDItools
        if grep -qE "\\[ERROR\\]|Killing job" \$LOG; then
            echo "REDItools a échoué pour $bam"
            exit 1
        fi

        # 🔍 Vérification output non vide (optionnel mais recommandé)
        if [ ! -s "\$OUT" ]; then
            echo "Output vide → erreur probable"
            exit 1
        fi
        """
}
