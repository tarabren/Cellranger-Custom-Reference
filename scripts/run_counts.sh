#!/bin/bash

###### PLUG: Define Variables #####
PROJECT_ID="_custom"  
PROJECT_DIRECTORY=""   
REFERENCE_DIRECTORY=""  
FASTQS=("")  # Array of FASTQ directories
BAM="true" 
##### NO PLUGS AFTER THIS LINE #####



# Initialize the master log file
MASTER_LOG="$PROJECT_DIRECTORY/$PROJECT_ID.log"  # Real-time log for all samples
echo "$(date '+%Y-%m-%d %H:%M:%S'): Starting job submission for project $PROJECT_ID" > "$MASTER_LOG"

# Ensure project, output, and error directories exist
mkdir -p "$PROJECT_DIRECTORY/outs"
mkdir -p "$PROJECT_DIRECTORY/err"

# Function to extract sample name from FASTQ file
extract_sample_name() {
    local file_name=$1
    # Extract everything before the first "_S[NUM]_" (Illumina naming convention)
    echo "$file_name" | sed -E 's/(_S[0-9]+_L[0-9]+_R[0-9]+_001\.fastq\.gz)$//'
}

# Identify all sample names from the FASTQ files in the FASTQ directories
SAMPLES=()
for fastq_dir in "${FASTQS[@]}"; do
    for fastq_file in "$fastq_dir"/*.fastq.gz; do
        if [[ -f "$fastq_file" && ! "$(basename "$fastq_file")" =~ ^Undetermined ]]; then
            sample_name=$(extract_sample_name "$(basename "$fastq_file")")
            # Avoid duplicate sample names
            if [[ ! " ${SAMPLES[@]} " =~ " ${sample_name} " ]]; then
                SAMPLES+=("$sample_name")
            fi
        fi
    done
done

# Check if there are no sample names found
if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "No FASTQ files found in the provided directories. Exiting." | tee -a "$MASTER_LOG"
    exit 1
fi

# Array to track job IDs
JOB_IDS=()

# Loop over each sample and check if both R1 and R2 files exist
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_LOG="$PROJECT_DIRECTORY/$SAMPLE.log"  # Log file for this specific sample

    # Collect R1 and R2 files across all FASTQ directories
    R1_FILES=()
    R2_FILES=()

    for fastq_dir in "${FASTQS[@]}"; do
        # Append matching R1 and R2 files from the current directory
        R1_FILES+=($(ls "$fastq_dir"/*"${SAMPLE}"*"_R1_"*.fastq.gz 2>/dev/null))
        R2_FILES+=($(ls "$fastq_dir"/*"${SAMPLE}"*"_R2_"*.fastq.gz 2>/dev/null))
    done

    # Start writing to the sample-specific log file
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Searching for FASTQs for sample $SAMPLE..." > "$SAMPLE_LOG"

    # Check if both R1 and R2 files were found and log them
    if [[ ${#R1_FILES[@]} -gt 0 ]] && [[ ${#R2_FILES[@]} -gt 0 ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Found the following R1 files for $SAMPLE:" >> "$SAMPLE_LOG"
        for r1 in "${R1_FILES[@]}"; do
            echo "  $r1" >> "$SAMPLE_LOG"
        done

        echo "$(date '+%Y-%m-%d %H:%M:%S'): Found the following R2 files for $SAMPLE:" >> "$SAMPLE_LOG"
        for r2 in "${R2_FILES[@]}"; do
            echo "  $r2" >> "$SAMPLE_LOG"
        done

        # Log the result to the master log
        echo "$(date '+%Y-%m-%d %H:%M:%S'): FASTQs found for $SAMPLE. See $SAMPLE_LOG for details." >> "$MASTER_LOG"
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Missing R1 or R2 files for sample $SAMPLE." >> "$SAMPLE_LOG"
        echo "$(date '+%Y-%m-%d %H:%M:%S'): Skipping $SAMPLE due to missing FASTQs." >> "$MASTER_LOG"
        continue
    fi

    # Create a temporary Slurm script for each sample
    SLURM_SCRIPT="${PROJECT_DIRECTORY}/submit_${SAMPLE}.sh"

    # Write the Slurm submission script for this sample
    cat > "$SLURM_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --nodes=1                  # Request 1 node for each job
#SBATCH --ntasks=1                 # 1 task per job
#SBATCH --cpus-per-task=16         # 16 CPUs per task
#SBATCH --mem=256g                 # Memory per node
#SBATCH --time=1-0:00              # Time limit
#SBATCH --job-name=${SAMPLE}_$PROJECT_ID  # Job name for each sample
#SBATCH --output=$PROJECT_DIRECTORY/outs/out_count_${SAMPLE}.txt  # Standard output file
#SBATCH --error=$PROJECT_DIRECTORY/err/err_count_${SAMPLE}.txt    # Standard error file

# Load cellranger module or other necessary environment setup
module load cellranger

# Run cellranger count for this sample
cellranger count --id="$SAMPLE" \
  --fastqs="${FASTQS[0]}" \
  --sample="$SAMPLE" \
  --transcriptome="$REFERENCE_DIRECTORY" \
  --create-bam="$BAM" \
  > "${PROJECT_DIRECTORY}/outs/out_count_${SAMPLE}.txt" \
  2> "${PROJECT_DIRECTORY}/err/err_count_${SAMPLE}.txt"

# Check if cellranger count succeeded
if [ \$? -ne 0 ]; then
  echo "cellranger count failed for sample $SAMPLE" | mail -s "Job Failure: $SAMPLE" "$NOTIFICATION_EMAIL"
  exit 1
fi

# Create success flag and send success notification
touch "${PROJECT_DIRECTORY}/success_${SAMPLE}.flag"
echo "cellranger count completed successfully for sample $SAMPLE" | mail -s "Job Success: $SAMPLE" "$NOTIFICATION_EMAIL"
EOF

    # Submit the script as a separate Slurm job and capture the job ID
    JOB_ID=$(sbatch "$SLURM_SCRIPT" | awk '{print $4}')  # Capture just the job ID
    JOB_IDS+=("$JOB_ID")  # Add the job ID to the array

    # Log the job submission time and job ID to the MASTER_LOG
    echo "$(date '+%Y-%m-%d %H:%M:%S'): Job for sample $SAMPLE submitted. Job ID: $JOB_ID" >> "$MASTER_LOG"
done

# Wait for all jobs to finish before submitting the aggregation job
echo "$(date '+%Y-%m-%d %H:%M:%S'): Waiting for all sample jobs to complete..." >> "$MASTER_LOG"
ALL_JOBS_COMPLETED=false
while [ "$ALL_JOBS_COMPLETED" == false ]; do
    ALL_JOBS_COMPLETED=true  # Assume all jobs are completed
    
    for JOB_ID in "${JOB_IDS[@]}"; do
        # Check if the job is still running in the queue
        JOB_STATUS=$(squeue -j "$JOB_ID" -h -o %T)
        if [[ "$JOB_STATUS" == "RUNNING" || "$JOB_STATUS" == "PENDING" ]]; then
            ALL_JOBS_COMPLETED=false
            break  # Exit the loop if any job is still running
        fi
    done
    
    if [ "$ALL_JOBS_COMPLETED" == false ]; then
        # Wait for a while before checking again (e.g., 10 minutes)
        sleep 600
    fi
done

# Once all count jobs are finished, submit the aggregation job
echo "$(date '+%Y-%m-%d %H:%M:%S'): All sample jobs completed. Submitting the aggregation job." >> "$MASTER_LOG"

# Create a separate job to handle aggregation after all count jobs finish
AGGREGATION_SCRIPT="${PROJECT_DIRECTORY}/submit_aggregation.sh"

cat > "$AGGREGATION_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128g
#SBATCH --time=1-0:00
#SBATCH --job-name=aggregation_$PROJECT_ID
#SBATCH --output=$PROJECT_DIRECTORY/outs/out_aggregation.txt
#SBATCH --error=$PROJECT_DIRECTORY/err/err_aggregation.txt

# Wait until all cellranger count jobs have completed successfully
MAX_ATTEMPTS=1440  # Max attempts (e.g., wait for up to 24 hours)
ATTEMPT=0
while true; do
    if [[ $(ls "$PROJECT_DIRECTORY"/success_*.flag 2>/dev/null | wc -l) -eq ${#SAMPLES[@]} ]]; then
        # Create aggr.csv
        echo "sample_id,molecule_h5" > "$PROJECT_DIRECTORY/aggr.csv"
        for SAMPLE in "${SAMPLES[@]}"; do
            echo "$SAMPLE,$PROJECT_DIRECTORY/$SAMPLE/outs/molecule_info.h5" >> "$PROJECT_DIRECTORY/aggr.csv"
        done

        echo "aggr.csv has been created in $PROJECT_DIRECTORY."

        # Run the cellranger aggr command
        if ! cellranger aggr --id=aggregated_analysis --csv="$PROJECT_DIRECTORY/aggr.csv" --normalize=none; then
            echo "cellranger aggr failed." | mail -s "Aggregation Job Failure" "$NOTIFICATION_EMAIL"
            exit 1
        fi

        # Send success notification
        echo "cellranger aggr completed successfully." | mail -s "Aggregation Job Success" "$NOTIFICATION_EMAIL"
        break
    fi
    ATTEMPT=$((ATTEMPT + 1))
    if [[ $ATTEMPT -ge $MAX_ATTEMPTS ]]; then
        echo "Not all jobs have completed within the expected time. Exiting." | mail -s "Aggregation Job Timeout" "$NOTIFICATION_EMAIL"
        exit 1
    fi
    sleep 600  # Check every 10 minutes
done
EOF

# Submit the aggregation job
sbatch "$AGGREGATION_SCRIPT"

echo "$(date '+%Y-%m-%d %H:%M:%S'): Job submission completed for project $PROJECT_ID" >> "$MASTER_LOG"