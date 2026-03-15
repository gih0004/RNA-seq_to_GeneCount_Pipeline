# RNA-seq_to_GeneCount_Pipeline
Automated RNA-seq pipeline designed for paired-sequencing for Quality Control Assesment (FASTQC), Adapter Trimming(Fastp),  Read Alignment(HISAT2), and Features[genes/cDNA/mRNA] Quantification (featurecounts) producing gene count matrices suitable for DESeq2 and other downstream bioinformatic analysis. 

The pipeline performs:    
1. Raw read quality control  
2. Adapter trimming and filtering  
3. Alignment to a reference genome  
4. Alignment quality summaries  
5. Gene read counting using featureCounts  
6. Generation of a gene count matrix suitable for differential expression analysis
 
The pipeline is designed to run on High Performance Computing (HPC) systems using SLURM.

## Required Input Files

Users must provide:
| File                     | Description                            | File Extension  |
| ------------------------ | -------------------------------------- |-----------------|
| FASTQ files              | Paired-end RNA-seq reads               | .fastq/fastq.gz |
| Reference genome FASTA   | Used for read alignment                |       .fa/.fna  |          
| Reference annotation GTF | Used for featureCounts gene assignment |   .gff3 / .gtf  |



## Software Dependencies
The pipeline requires the following bioinformatics tools:  
| Software                | Purpose                             |
|-------------------------|-------------------------------------|
| FastQC                  | Raw read quality assessment         |
| fastp                   | Adapter trimming and read filtering |
| HISAT2                  | RNA-seq read alignment              |
| SAMtools                | SAM/BAM file processing             |
| Subread (featureCounts) | Gene expression quantification      |
| MultiQC                 | Aggregated quality control reporting|



## Directory Setup
To run the pipeline, create a main directory containing:   
A. Reference genome FASTA file  
B. Reference annotation GTF/GFF3 file  
C. The pipeline script  
D. Pipeline Submission Script  

Make Subdirectory titled <ins> RNA_Reads</ins>  containing:  
E. Raw FASTQ sequencing files 

### FASTQ File Naming Requirements
The pipeline expects **paired-end FASTQ files** following a specific naming convention so that forward and reverse reads can be detected automatically.
Files must follow this format:

|         Example             |  Suffix   |      Description       |
|-----------------------------|-----------|------------------------|
|   sample_name_1.fastq       |    `_1`   | Forward reads (Read 1) |
|   sample_name_2.fastq       |    `_2`   | Reverse reads (Read 2) |

> Compressed FASTQ files are also supported: sample_name_1.fastq.gz , sample_name_2.fastq.gz


## Directory diagram  
```ruby
>Project_Directory/
```
* reference_genome.fa              
* annotation.gtf                
* GHC_RNA_pipeline.sh                 
* GHC_Pipeline_submission.sh  ]
```ruby      
>Project_Directory/RNA_Reads          # Subdirectory containing raw sequencing reads
```
>[!Warning]
>RNA_Reads MUST be the name of the directory

```ruby
>Project_Directory/Results            #This is produced by the pipeline 
```

-----------------------------------------------------

## Running this Pipeline
There are two ways to run the pipeline.
1. Use the Submission Helper Script (Recommended) - GHC_Pipeline_Submission.sh
```ruby
$ bash 'GHC_Pipeline_Submission.sh'
```
The script will:
1. Automatically detect reference genome FASTA files in the directory
2. Automatically detect GTF/GFF annotation files
3. Prompt the user to enter a species name
4. Display the most common feature types present in the annotation file
5. Ask which feature type should be counted
6. Submit the pipeline to SLURM


2. Manual Submission 
Advanced users may submit the pipeline directly using sbatch.
```ruby
$ sbatch GHC_RNA_Pipeline.sh <species> <reference_fasta> <reference_gtf> <feature_type> <gene_attribute>

#Example:
#$ sbatch GHC_RNA_Pipeline.sh arabidopsis TAIR10.fa TAIR10.gtf mRNA
```

| Argument        | Description                                                                                 |  
| --------------- | --------------------------------------------------------------------------------------------|
| species         | Name used to build the HISAT2 index                                                         |
| reference_fasta | Reference genome FASTA file                                                                 |
| reference_gtf   | Annotation file used by featureCounts                                                       |
| feature_type    | Feature type counted by featureCounts (e.g. `gene`, `mRNA`, `CDS`, `exon`)                  |
| gene_attribute  |Attribute from the GTF/GFF file used by featureCounts to group reads into genes (`-g` option)| 

> If gene attribute is omitted, the pipeline automatically selects a compatible attribute (e.g., `ID` for `gene`, `Parent` for transcript features).
-----------------------------------------------------
### Pipeline Logging

Each pipeline run automatically generates a timestamped log file that records pipeline progress and errors.  
Example log file: GHC_RNA_Pipeline_20260125_143210.log  

The log file contains:
1. Pipeline start time and arguments used
2. Progress updates for each pipeline step
3. Alignment and processing status
4. Error messages if a step fails  
This file is written to the /Project_Directory and can be monitored during execution to track pipeline progress.

-----------------------------------------------------



### HPC Configuration Example

These are the configurations I usewhen running this script on my HPC
```ruby 
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=6-0
```
Adjust these parameters depending on dataset size and cluster configuration.

-----------------------------------------------------

# Pipeline Workflow

### Step 1A — Quality Control of Raw Reads

The pipeline first evaluates raw sequencing read quality using FastQC.
```ruby
fastqc "$RAW_READS_DIR"/*.fastq* -o "$FASTQC_DIR"
```
>Results are stored in:
>```ruby
> Project_Directory/Results/FASTQC
>```

For information on FASTQC results, please see:  
[Results Analysis Video](https://www.youtube.com/watch?v=bz93ReOv87Y)  
[Results Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

-----------------------------------------------------

### Step 1B - Adapter Removal and Read Filtering
Reads are filtered using fastp to remove adapters and low-quality sequences.
```ruby
fastp \
-i sample_1.fastq \
-I sample_2.fastq \
-o sample.filtered.1.fq \
-O sample.filtered.2.fq
```

fastp also produces quality reports:
```ruby
sample_fastp.html
sample_fastp.json
```

Filtered reads are stored in:
```ruby
>Project_Directory/Results/FASTQC_filtered/
```


----------------------------------------
### Step 1C — Quality Control of Filtered Reads

FastQC is run again on filtered reads to confirm quality improvement.

```ruby
fastqc "$FASTQC_FILTERED_DIR"/*.filtered.1.fq -o "$FASTQC_FILTERED_DIR"
fastqc "$FASTQC_FILTERED_DIR"/*.filtered.2.fq -o "$FASTQC_FILTERED_DIR"
```

Individual FASTQC Html results for FASTQC are moved into specific directories
```ruby
mv "$FASTQC_DIR"/*.html "$HTML_DIR"
mv "$FASTQC_FILTERED_DIR"/*.html "$HTML_FILTERED_DIR"
```
>Namely:
> For unfiltered samples QC:
>```ruby
>>Project_Directory/Results/FASTQC_html
>```
>For filtered samples QC:
>```ruby
>>Project_Directory/Results/FASTQC_Filtered_html
>```


There is also an aggregated QC report for easy viewing which is produced in the main directory
```ruby
multiqc $FASTQC_FILTERED_DIR/. -o "$FASTQC_FILTERED_DIR"
>Project_Directory/Results/FASTQC_filtered/multiqc_report.html
```
---------------------------------------------------

### STEP 2: HISAT2 indexing

If a HISAT2 index does not exist, the pipeline builds one

```ruby 
hisat2-build "$reference_fa" "$BASE_DIR/$species"
log "Finished building HISAT2 index"
```
>[!TIP]
>Note that most steps are being logged. For progress updates when the pipeline is running, you can view the log file, its purpose is to document where within the pipeline is the HPC at that time.

The indexed genome files generated will be in the main directory and look like:
```ruby
species.1.ht2
species.2.ht2
...
```
------------------------------------------------

### Step 3: HISAT2 Read Alignment
After creating indices from genome, we can then run alignment of the samples against the indices recently created:  

```ruby   
    hisat2 -p $THREADS -q -x "$BASE_DIR/$species" \
        -1 "$fq1" -2 "$fq2" \
        -S "$SAM_DIR/${sample_name}_aligned.sam" \
        --summary-file "$ALIGN_SUMMARY_DIR/${sample_name}_summary.txt"

```
> When using for big datasets, threads can be changed under "Thread Configuration" Section in GHC_RNA_Pipeline.sh code 

Alignment Summaries are written to: 
```ruby
> Project_Directory/Results/Sample_Alignments_Summary
```

For more information on HISAT2 read alignment options and alignment summaries please see:  
[HISAT2 Documentation](https://daehwankimlab.github.io/hisat2/manual/)

--------------------------------------------------
  
### Step 4: Alignment Summary Table
The pipeline generates a summary table combining relevant alignment statistics.

|    Sample   |   Total Reads  | Overall Alignment Score (%) | Aligned 0 times concordantly |
| ----------- | ---------------|-----------------------------|------------------------------|
|   Sample 1  |    11758107    |          99.06              |         226868 (1.93%)       |

>To view Alignment Summary as the pipeline runs, use `Alignment_Summary_Table.txt`.   
>This tab delimited table is downloadable and produced in the main project directory
```ruby
>Project_Directory/Alignment_Summary_Table
```

-----------------------------------------------------

## Step 5: Convert SAM Files to BAM Files
To do anything meaningful with alignment data you must swith from SAM to its binary counterpart BAM file.  
This binary format is much easier for computer programs to work with.

SAM files are converted into sorted BAM files using samtools.

```ruby
samtools sort sample.sam -o sample.bam
samtools index sample.bam
```
>[!Note]
> SAM files are removed after conversion to reduce disk usage.
>```ruby
>rm -r $SAM_DIR
>```
Final BAM files are stored in:
```ruby
>Project_Directory/Results/BAM_Files
```
----------------------------------------------------

## STEP 6: Run featurecounts to generate gene counts 
Gene expression counts are generated using featureCounts, a program from the Subread package that assigns aligned reads to genomic features defined in the annotation file.

```ruby 
    featureCounts -p -O -T $THREADS -g "$gene_attribute" -f -F GTF -t "$feature_type" -M --countReadPairs \
        -a "$reference_gtf" \
        -o "$FEATURECOUNTS_DIR/${sample_name}.txt"  $bam
```
>[!Important]
>The feature type ($feature_type) is selected when running the submission script.  
>This value corresponds directly to the -t option in featureCounts, which specifies which genomic feature type should be counted from the annotation file.

>Example command during pipeline submission:  
>Format: GHC_RNA_Pipeline.sh <species> <reference_fasta> <reference_gtf> <feature_type> <gene_attribute>
>```ruby
>$ sbatch GHC_RNA_Pipeline.sh Arabidopsis TAIR10.fa TAIR10.gtf mRNA 
>```
>In this example:  
>```ruby
>-t mRNA
>```
>tells featureCounts to count reads overlapping mRNA features in the GTF file.

##### Feature Types in a GTF/GFF File
| Feature Type | Description              |
| ------------ | ------------------------ |
| gene         | Entire gene region       |
| mRNA         | Messenger RNA transcript |
| exon         | Individual exons         |
| CDS          | Coding sequence          |
| transcript   | Transcript features      |

The submission helper script automatically displays the most common feature types found in the annotation file  
This allows the user to select the correct feature type for their analysis.


#### Pipeline featureCounts Options:
| Option             | Description                                              |
| ------------------ | -------------------------------------------------------- |
| `-p`               | Indicates paired-end sequencing data                     |
| `-O`               | Allows reads overlapping multiple features to be counted |
| `-T`               | Number of threads used                                   |
| `-g`               | Attribute used to group features (gene identifier)       |
| `-f`               | Counts reads at the feature level                        |
| `-F GTF`           | Specifies that the annotation format is GTF              |
| `-t feature_type`  | Feature type to count (selected during submission)       |
| `-M`               | Counts multi-mapping reads                               |
| `--countReadPairs` | Counts fragments instead of individual reads             |

Complete documentation for featureCounts can be found in its Subread Manual:  
[featureCounts Documentation](https://subread.sourceforge.net/SubreadUsersGuide.pdf)

To prevent mismatches between -t (feature type) and -g (gene identifier), the pipeline automatically selects a compatible attribute based on the chosen feature type.  
Typical pairings used by the pipeline:

| Feature Type (`-t`) | Identifier (`-g`) |
| ------------------- | ----------------- |
| gene                | ID                |
| mRNA / transcript   | Parent            |
| exon                | Parent            |
| CDS                 | Parent            |


This ensures reads are grouped correctly at the gene level rather than incorrectly counting individual transcript features.  
The selected attribute is printed in the pipeline log file.  

Individual Sample Count files are stored in:
```ruby
>Project_Directory/Results/Featurecounts
```

------------------------------------------------------------
## Step 7: Gene Count Matrix Generation

After featureCounts has generated individual count files for each sample, the pipeline combines these results into a single gene count matrix.  
This matrix is the primary output used for downstream RNA-seq analysis tools such as DESeq2, edgeR, or limma.
|    Gene   |   Sample A   | Sample B | Sample C |
| --------- | -------------|----------|----------|
|   Gene 1  |     26       |  136     |    7     |
|   Gene 2  |     8        |  23      |    91    |
|   Gene 3  |     6        |  254     |    2     |

The final matrix is written to:
```ruby
>Project_Directory/Gene_Count_Matrix.txt
```


-------------------------------------------------------------
### Expected Outputs
After successful completion the pipeline generates:
| File                          | Description                              |
| ----------------------------- | ---------------------------------------- |
| `Alignment_Summary_Table`     | Alignment statistics for all samples     |
| `Results/BAM_Files/*.bam`     | Sorted BAM alignment files               |
| `Results/Featurecounts/*.txt` | Raw featureCounts output for each sample |
| `Gene_Count_Matrix.txt`       | Combined gene expression matrix          |
| `RNA_Pipeline_TIMESTAMP.log`  | Pipeline execution log                   |



