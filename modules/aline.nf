/*
Adapted from
from https://github.com/mahesh-panchal/nf-cascade
Auto-detects HPC/local environment and adapts execution strategy
*/
process AliNe {
    tag "$pipeline_name"
    label 'aline'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val pipeline_name     // String
        val profile           // String
        val config
        val reads
        val genome
        val read_type
        val aligner
        val library_type
        val annotation
        val cache_dir          // String
        val use_slurm          // Boolean - whether parent is using slurm

    when:
        task.ext.when == null || task.ext.when

    script:
    def nxf_cmd = "nextflow run ${pipeline_name} ${profile} ${config} --reads ${reads} --reference ${genome} ${read_type} ${aligner} ${library_type} --annotation ${annotation} --data_type rna --outdir ${task.workDir}/AliNe"
    """
    # Create cache directory
    mkdir -p "${cache_dir}"
    cd "${cache_dir}"
    
    # Save command for reference/debugging
    echo "${nxf_cmd}" > ${task.workDir}/nf-cmd.sh
    
    # Detect execution environment and run AliNe accordingly
    if ${use_slurm} && command -v sbatch >/dev/null 2>&1; then
        echo "[AliNe] Detected HPC environment - submitting AliNe as separate SLURM job"
        
        # Create sbatch script for AliNe
        cat > ${task.workDir}/aline_job.sh <<'SBATCH_EOF'
#!/bin/bash
#SBATCH --job-name=AliNe_pipeline
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=${task.workDir}/aline_%j.out
#SBATCH --error=${task.workDir}/aline_%j.err

set -euo pipefail

cd "${cache_dir}"

echo "Starting AliNe pipeline at \$(date)"
${nxf_cmd}
echo "AliNe pipeline completed at \$(date)"
SBATCH_EOF
        
        # Submit job and capture job ID
        JOB_ID=\$(sbatch --parsable ${task.workDir}/aline_job.sh)
        echo "[AliNe] Submitted SLURM job: \$JOB_ID"
        echo \$JOB_ID > ${task.workDir}/aline_job_id.txt
        
        # Wait for job completion
        echo "[AliNe] Waiting for job \$JOB_ID to complete..."
        while squeue -j \$JOB_ID 2>/dev/null | grep -q \$JOB_ID; do
            sleep 30
        done
        
        # Check job exit status
        JOB_STATE=\$(sacct -j \$JOB_ID --format=State --noheader | head -1 | tr -d ' ')
        echo "[AliNe] Job \$JOB_ID finished with state: \$JOB_STATE"
        
        if [[ "\$JOB_STATE" != "COMPLETED" ]]; then
            echo "[AliNe] ERROR: Job failed with state \$JOB_STATE" >&2
            cat ${task.workDir}/aline_*.err >&2 || true
            exit 1
        fi
        
        # Copy log for reference
        if [ -f .nextflow.log ]; then
            cp .nextflow.log ${task.workDir}/nextflow.log
        fi
        
        echo "[AliNe] Pipeline completed successfully via SLURM"
    else
        echo "[AliNe] Detected local/standard environment - running AliNe directly"
        
        # Run nextflow command directly
        ${nxf_cmd}
        
        # Copy log for reference
        if [ -f .nextflow.log ]; then
            cp .nextflow.log ${task.workDir}/nextflow.log
        fi
        
        echo "[AliNe] Pipeline completed successfully (direct execution)"
    fi
    """

    output:
        path "AliNe"  , emit: output
}