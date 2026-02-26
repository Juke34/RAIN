/*
Adapted from
from https://github.com/mahesh-panchal/nf-cascade
Auto-detects HPC/local environment and adapts execution strategy
*/
process AliNe {
    tag "$pipeline_name"
    label 'aline'
    publishDir "${params.outdir}", mode: 'copy'
    
    maxRetries 0  // Override global retry config - do not retry this process

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
        def nxf_cmd = "nextflow run ${pipeline_name} ${profile} ${config} --reads ${reads} --reference ${genome} ${read_type} ${aligner} ${library_type} --annotation ${annotation} --data_type rna --outdir \$WORK_DIR/AliNe"
        """
        echo "[AliNe] Process started at \$(date '+%Y-%m-%d %H:%M:%S')"
        
        # Save absolute work directory before changing context
        WORK_DIR=\$(pwd)

        # Create cache directory for resume AliNe run made from different working directory
        mkdir -p "${cache_dir}"
        cd "${cache_dir}"

        # Save command for reference/debugging
        echo "${nxf_cmd}" > \$WORK_DIR/nf-cmd.sh

        # Detect execution environment and run AliNe accordingly
        if ${use_slurm} && command -v sbatch >/dev/null 2>&1; then
            echo "[AliNe] Detected HPC environment - submitting AliNe as separate SLURM job"
            
            # Create sbatch script for AliNe
            cat > \$WORK_DIR/aline_job.sh <<SBATCH_EOF
        #!/bin/bash
        #SBATCH --job-name=rain_AliNe_pipeline
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=4G
        #SBATCH --time=2-00:00:00
        #SBATCH --output=\$WORK_DIR/aline_%j.out
        #SBATCH --error=\$WORK_DIR/aline_%j.err

        set -euo pipefail

        cd "${cache_dir}"

        echo "Starting AliNe pipeline at \$(date)"
        ${nxf_cmd}
        echo "AliNe pipeline completed at \$(date)"
        SBATCH_EOF
            
            # Submit job and capture job ID
            JOB_ID=\$(sbatch --parsable \$WORK_DIR/aline_job.sh)
            echo "[AliNe] Submitted SLURM job: \$JOB_ID at \$(date '+%Y-%m-%d %H:%M:%S')" 
            echo \$JOB_ID > \$WORK_DIR/aline_job_id.txt
            
            # Wait for job to appear in scheduler queue
            echo "[AliNe] Waiting for job to appear in scheduler queue..."
            sleep 5
            
            RETRY=0
            while [ \$RETRY -lt 12 ]; do
                if squeue -j \$JOB_ID 2>/dev/null | grep -q \$JOB_ID; then
                    echo "[AliNe] Job \$JOB_ID is now visible in queue"
                    break
                fi
                echo "[AliNe] Job not yet visible, waiting... (attempt \$((RETRY+1))/12)"
                sleep 5
                RETRY=\$((RETRY+1))
            done
            
            # Wait for job completion
            echo "[AliNe] Waiting for job \$JOB_ID to complete..."
            while squeue -j \$JOB_ID 2>/dev/null | grep -q \$JOB_ID; do
                sleep 30
            done
            
            # Check job exit status
            JOB_STATE=\$(sacct -j \$JOB_ID --format=State --noheader | head -1 | tr -d ' ')
            echo "[AliNe] Job \$JOB_ID finished with state: \$JOB_STATE at \$(date '+%Y-%m-%d %H:%M:%S')"
            
            if [[ "\$JOB_STATE" != "COMPLETED" ]]; then
                echo "[AliNe] ERROR: Job failed with state \$JOB_STATE at \$(date '+%Y-%m-%d %H:%M:%S')" >&2
                cat \$WORK_DIR/aline_*.err >&2 || true
                exit 1
            fi
            
            # Copy log for reference
            if [ -f .nextflow.log ]; then
                cp .nextflow.log \$WORK_DIR/nextflow.log
            fi
            
            echo "[AliNe] Pipeline completed successfully via SLURM at \$(date '+%Y-%m-%d %H:%M:%S')"
        else
            echo "[AliNe] Detected local/standard environment - running AliNe directly"
            
            # Run nextflow command directly
            ${nxf_cmd} || {
                echo "[AliNe] ERROR: Pipeline failed at \$(date '+%Y-%m-%d %H:%M:%S')" >&2
                exit 1
            }
            
            # Copy log for reference
            if [ -f .nextflow.log ]; then
                cp .nextflow.log \$WORK_DIR/nextflow.log
            fi
            
            echo "[AliNe] Pipeline completed successfully (direct execution) at \$(date '+%Y-%m-%d %H:%M:%S')"
        fi
        
        echo "[AliNe] Process finished at \$(date '+%Y-%m-%d %H:%M:%S')"
"""

    output:
        path "AliNe"  , emit: output
}