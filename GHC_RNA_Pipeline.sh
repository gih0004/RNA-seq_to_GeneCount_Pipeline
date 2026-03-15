#!/bin/bash

# ==============================
# SLURM JOB CONFIGURATION
# ==============================

#SBATCH --nodes=1                        # Number of compute nodes requested
#SBATCH --mem=64G                        # Total RAM requested
#SBATCH --time=6-0                       # Max runtime (6 days)


# ==============================
# SHELL SAFETY
# ==============================

set -euo pipefail


# ==============================
# LOGGING CONFIGURATION
# ==============================

BASE_DIR=$(pwd)

LOGFILE="$BASE_DIR/RNA_Pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}


# ==============================
# ERROR HANDLING + DEFAULTS
# ==============================

trap 'log "Pipeline crashed at line $LINENO"' ERR

log "Pipeline started"
log "Default threads used for hisat2 and featureCounts commands are set at 1"


# ==============================
# SOFTWARE MODULE LOADS
# ==============================
# Loads required bioinformatics tools from the HPC module system
log "Loading software modules"

module load FastQC/0.12.1-Java-11        # Read quality assessment
module load HISAT2/2.2.1-gompi-2023a     # RNA-seq aligner
module load StringTie/2.2.3-GCC-12.3.0   # Transcript assembly (not used here yet)
module load SAMtools/1.18-GCC-12.3.0     # SAM/BAM processing
module load fastp/0.23.4-GCC-12.3.0      # Read filtering and trimming
module load Subread/2.0.6-GCC-12.3.0     # Provides featureCounts
module load MultiQC/1.22.3-foss-2023b    # Aggregate QC reports

log "Modules loaded successfully"

# ==============================
# PIPELINE CONFIGURATION
# ==============================

# Input directory containing raw FASTQ reads
RAW_READS_DIR="$BASE_DIR/RNA_Reads" #Directory should be created by user

#Results Directory
RESULTS_DIR="$BASE_DIR/Results"

# QC output directories
FASTQC_DIR="$RESULTS_DIR/FASTQC"  #FASTQC files of unfiltered samples
FASTQC_FILTERED_DIR="$RESULTS_DIR/FASTQC_filtered" #FASTQC files of filtered samples

# HTML report storage directories
HTML_DIR="$RESULTS_DIR/FASTQC_html"    #Downloadable html reports for samples
HTML_FILTERED_DIR="$RESULTS_DIR/FASTQC_Filtered_html" #Downloadable html reports for filtered samples

# Alignment summary directory
ALIGN_SUMMARY_DIR="$RESULTS_DIR/Sample_Alignments_Summary" #Individual sample summary tables

# featureCounts output directory
FEATURECOUNTS_DIR="$RESULTS_DIR/Featurecounts"

# Alignment file storage
SAM_DIR="$RESULTS_DIR/SAM_Files"
BAM_DIR="$RESULTS_DIR/BAM_Files"

# ==============================
# INPUT PARAMETERS
# ==============================
# Assign CLI inputs to variables
species=$1
reference_fa=$2
reference_gtf=$3
feature_type=$4
gene_attribute=${5:-}


#==============================
# INPUT ARGUMENT VALIDATION
# ==============================
# Expected arguments:
#   1. Species name (used for index naming)
#   2. Reference genome FASTA
#   3. Reference annotation GTF
#   4. Feature type to count (e.g., CDS, exon, gene)
#   5. Feature ID to use for counts (e.g.  Parent, ID, Transcript)

#Validate arguments
if [ "$#" -lt 4 ]; then
  log "Usage: $0 <species> <reference_fasta> <reference_gtf> <feature_type>"
  exit 1
fi

# Validate input files exist
if [ ! -f "$reference_fa" ]; then
    log "ERROR: Reference genome not found: $reference_fa"
    exit 1
fi

if [ ! -f "$reference_gtf" ]; then
    log "ERROR: Reference GTF not found: $reference_gtf"
    exit 1
fi

log "Scanning for FASTQ files in $RAW_READS_DIR"
if [ ! -d "$RAW_READS_DIR" ]; then
    log "ERROR: Raw reads directory not found: $RAW_READS_DIR"
    exit 1
fi

# Auto-detect gene attribute if not provided
if [[ -z "$gene_attribute" ]]; then
    if [[ "$feature_type" == "gene" ]]; then
        gene_attribute="ID"
    else
        gene_attribute="Parent"
    fi
fi


# ==============================
# INPUT CONFIGURATION SUMMARY
# ==============================
# Print configuration summary
log "RNA-seq Pipeline"
log "Run directory: $BASE_DIR"
log "Species: $species"
log "Reference Genome: $reference_fa"
log "Reference GTF: $reference_gtf"
log "Feature type to count: $feature_type"
log "featureCounts ID to use:$gene_attribute"

# ==============================
# DIRECTORY CREATION
# ==============================
# Creates all required output directories if they do not already exist

mkdir -p "$FASTQC_DIR" \
         "$FASTQC_FILTERED_DIR" \
         "$HTML_DIR" \
         "$HTML_FILTERED_DIR" \
         "$ALIGN_SUMMARY_DIR" \
         "$FEATURECOUNTS_DIR" \
         "$SAM_DIR" \
         "$BAM_DIR"


# ==============================
# STEP 1A — FASTQC ON RAW READS
# ==============================
log "Scanning for FASTQ files in $RAW_READS_DIR"

sample_count=$(ls "$RAW_READS_DIR"/*1.fastq* 2>/dev/null | wc -l)

log "Detected $sample_count paired-end samples"

# Generates quality reports for raw FASTQ files
log "STEP 1A: Running FastQC on raw unfiltered reads "
fastqc "$RAW_READS_DIR"/*.fastq* -o "$FASTQC_DIR"


# ==============================
# STEP 1B — READ FILTERING (fastp)
# ==============================
# Adapter trimming, quality filtering, and report generation
log "STEP 1B: Starting read filtering with fastp"

for fq1 in "$RAW_READS_DIR"/*1.fastq*
do
    base="${fq1%1.fastq*}"                     # Remove read1 suffix to locate paired read2
    sample_name=$(basename "${fq1%.fastq*}")   # Extract sample name
    log "Filtering reads for sample: $sample_name"
    fastp \
        -i "$fq1" \
        -I ${base}2.fastq* \
        -o $FASTQC_FILTERED_DIR/${sample_name}.filtered.1.fq \
        -O $FASTQC_FILTERED_DIR/${sample_name}.filtered.2.fq \
        --detect_adapter_for_pe --qualified_quality_phred 20 \
        -h $FASTQC_FILTERED_DIR/${sample_name}_fastp.html \
        -j $FASTQC_FILTERED_DIR/${sample_name}_fastp.json
done


# ==============================
# STEP 1C — FASTQC ON FILTERED READS
# ==============================
# Quality control after trimming/filtering
log "STEP 1C: Running FastQC on filtered reads"

fastqc "$FASTQC_FILTERED_DIR"/*.filtered.1.fq -o "$FASTQC_FILTERED_DIR"
fastqc "$FASTQC_FILTERED_DIR"/*.filtered.2.fq -o "$FASTQC_FILTERED_DIR"

# Aggregates QC reports
multiqc $FASTQC_FILTERED_DIR/. -o "$FASTQC_FILTERED_DIR"


# Copy HTML reports into dedicated directories
mv "$FASTQC_DIR"/*.html "$HTML_DIR"
mv "$FASTQC_FILTERED_DIR"/*.html "$HTML_FILTERED_DIR"


# ==============================
# THREAD CONFIGURATION
# ==============================
# Default thread count used by HISAT2
THREADS=1


# ==============================
# STEP 2 — HISAT2 INDEXING
# ==============================
# Builds genome index if it does not already exist
log "STEP 2: Checking HISAT2 index"
if [ ! -f "$BASE_DIR/$species.1.ht2" ]; then
    log "Building HISAT2 index for $species"
    hisat2-build "$reference_fa" "$BASE_DIR/$species"
    log "Finished building HISAT2 index"
else
    log "Existing HISAT2 index detected"
fi

# ==============================
# STEP 3 — READ ALIGNMENT (HISAT2)
# ==============================
# Aligns paired-end reads to the reference genome
log "STEP 3: Starting HISAT2 alignment"
for fq1 in "$FASTQC_FILTERED_DIR"/*.filtered.1.fq
do
    sample_name=$(basename "${fq1%.filtered.1.fq}")    # Extract sample name
    fq2="${FASTQC_FILTERED_DIR}/${sample_name}.filtered.2.fq"   # Paired read
    
    log "Aligning sample: $sample_name" 
    
    hisat2 -p $THREADS -q -x "$BASE_DIR/$species" \
        -1 "$fq1" -2 "$fq2" \
        -S "$SAM_DIR/${sample_name}_aligned.sam" \
        --summary-file "$ALIGN_SUMMARY_DIR/${sample_name}_summary.txt"

    log "Finished hisat2 alignment $sample_name" 
done 


# ==============================
# STEP 4 — ALIGNMENT SUMMARY TABLE
# ==============================
# Extracts key metrics from HISAT2 summary outputs
log "STEP 4: Generating alignment summary table for all samples"


SUMMARY_FILE="$BASE_DIR/Alignment_Summary_Table.txt"

# Header
echo -e "File\tTotal Reads\tAlignment Score (%)\tPairs Aligned Concordantly 0 times" > "$SUMMARY_FILE"
log "Compiling alignment summary statistics"

for file in "$ALIGN_SUMMARY_DIR"/*_summary.txt
do

    sample_name=$(basename "${file%_summary.txt}")

    # Total reads from first line
    total_reads=$(head -1 "$file" | grep -o '^[0-9]\+')

    # Overall alignment percentage
    alignment_score=$(tail -n 1 "$file" | grep -oP '\d+\.\d+(?=%)')

    # Number of pairs failing concordant alignment
    pairs_aligned_0=$(grep "aligned concordantly 0 times" "$file" | head -n 1 | awk '{print $1" "$2}')

    echo -e "$sample_name\t$total_reads\t$alignment_score\t$pairs_aligned_0" >> "$SUMMARY_FILE"
done


# ==============================
# STEP 5 — SAM → SORTED BAM
# ==============================
# Converts SAM alignments into sorted BAM files
log "STEP 5: Converting SAM to sorted BAM"

for sam in "$SAM_DIR"/*_aligned.sam
do
    sample_name=$(basename "${sam%_aligned.sam}")
    log "Converting $sample_name SAM → BAM"
    samtools sort "$sam" -o "$BAM_DIR/${sample_name}.bam" -O bam
    samtools index "$BAM_DIR/${sample_name}.bam"
    

    log "Finished SAM to BAM conversion of $sample_name"
done

rm -r $SAM_DIR

# ==============================
# STEP 6 — FEATURECOUNTS
# ==============================
# Quantifies reads mapped to genomic features
log "STEP 6: Running featureCounts"
for bam in "$BAM_DIR"/*.bam
do
    sample_name=$(basename "${bam%.bam}")
    log "Running featureCounts on $sample_name"
    featureCounts -p -O -T $THREADS -g "$gene_attribute" -f -F GTF -t "$feature_type" -M --countReadPairs \
        -a "$reference_gtf" \
        -o "$FEATURECOUNTS_DIR/${sample_name}.txt"  $bam

    # Remove temporary files created by featureCounts
    rm -f *temp

    log "Finished featurecounts $sample_name processing" 
done

# ==============================
# STEP 7: Gene Count Matrix Table
# ==============================

log "STEP 7: Combining featureCounts outputs"

COUNTS_MATRIX="$BASE_DIR/Gene_Count_Matrix.txt"

first_file=$(ls "$FEATURECOUNTS_DIR"/*.txt | head -n 1)

# Extract gene IDs from first file
cut -f1 "$first_file" | grep -v '^#' > "$COUNTS_MATRIX"

#Sanity Check in case results are not found to avoid silent crashing
count_files=$(ls "$FEATURECOUNTS_DIR"/*.txt 2>/dev/null | wc -l)

if [ "$count_files" -eq 0 ]; then
    log "ERROR: No featureCounts output files found"
    exit 1
fi

# Add count columns from each sample
for file in "$FEATURECOUNTS_DIR"/*.txt
do
    sample=$(basename "$file" .txt)
    cut -f7 "$file" | grep -v '^#' > "$FEATURECOUNTS_DIR/${sample}.counts"
done

paste "$COUNTS_MATRIX" "$FEATURECOUNTS_DIR"/*.counts > temp_matrix
mv temp_matrix "$COUNTS_MATRIX"

log "Combining $(ls $FEATURECOUNTS_DIR/*.txt | wc -l) samples into gene count matrix"
log "Gene count matrix created: $COUNTS_MATRIX"
rm "$FEATURECOUNTS_DIR"/*.counts
log "RNA-seq pipeline completed successfully"
log "Results stored in $BASE_DIR"