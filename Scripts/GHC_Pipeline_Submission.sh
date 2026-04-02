#!/bin/bash

echo "=== RNA-seq Pipeline Submission ==="

# Auto-detect FASTA
default_fa=$(ls ./*.fa* 2>/dev/null | head -n 1)
if [[ -z "$default_fa" ]]; then
  read -p "No FASTA file found. Enter path to reference FASTA: " reference_fa
else
  read -p "Enter reference FASTA path [$default_fa]: " reference_fa
  reference_fa="${reference_fa:-$default_fa}"
fi

# Auto-detect GTF or GFF3
default_gtf=$(ls ./*.gtf* ./*.gff* 2>/dev/null | head -n 1)
if [[ -z "$default_gtf" ]]; then
  read -p "No GTF/GFF3 file found. Enter path to annotation: " reference_gtf
else
  read -p "Enter reference GTF/GFF3 path [$default_gtf]: " reference_gtf
  reference_gtf="${reference_gtf:-$default_gtf}"
fi

# Ask for species name
read -p "Enter species name (e.g., Arabidopsis_thaliana): " species

# Show available feature types
echo -e "\nAvailable feature types in $reference_gtf:"
cut -f3 "$reference_gtf" | grep -v '^#' | sort | uniq -c | sort -nr | head -n 10

# Prompt for feature type
read -rp "Enter feature type to count (e.g., exon, CDS, mRNA, gene): " feature_type

# Extract attribute keys from 9th column (GFF/GTF)
mapfile -t attributes < <(
    cut -f9 "$reference_gtf" \
    | tr ';' '\n' \
    | cut -d'=' -f1 \
    | cut -d' ' -f1 \
    | sort -u
)

# Auto-detect default
if [[ -z "$gene_attribute" ]]; then
    if [[ "$feature_type" == "gene" ]]; then
        default_attr="ID"
    elif printf '%s\n' "${attributes[@]}" | grep -qx "Parent"; then
        default_attr="Parent"
    else
        default_attr="ID"
    fi
fi

echo "Detected attributes: ${attributes[*]}"
echo "Suggested gene attribute: $default_attr"


# Validate
if ! printf '%s\n' "${attributes[@]}" | grep -qx "$gene_attribute"; then
    echo "ERROR: '$gene_attribute' not found in annotation attributes."
    exit 1
fi

echo "Using gene_attribute=$gene_attribute"

\
# Confirm
echo -e "\nSubmitting pipeline with:"
echo "Species        : $species"
echo "Reference FASTA: $reference_fa"
echo "Annotation GTF : $reference_gtf"
echo "Feature Type   : $feature_type"
echo "Gene Attribute : $gene_attribute"


# Final confirmation
read -p "Proceed with submission? [Y/n]: " confirm
confirm=${confirm:-Y}

if [[ "$confirm" =~ ^[Yy]$ ]]; then
    sbatch GHC_RNA_Pipeline.sh "$species" "$reference_fa" "$reference_gtf" "$feature_type" "$gene_attribute"
else
    echo "Submission cancelled."
fi
echo "Using gene_attribute=$gene_attribute"


# Confirm
echo -e "\nSubmitting pipeline with:"
echo "Species        : $species"
echo "Reference FASTA: $reference_fa"
echo "Annotation GTF : $reference_gtf"
echo "Feature Type   : $feature_type"
echo "Gene Attribute : $gene_attribute"


# Final confirmation
read -p "Proceed with submission? [Y/n]: " confirm
confirm=${confirm:-Y}

if [[ "$confirm" =~ ^[Yy]$ ]]; then
    sbatch GHC_RNA_Pipeline.sh "$species" "$reference_fa" "$reference_gtf" "$feature_type" "$gene_attribute"
else
    echo "Submission cancelled."
fi