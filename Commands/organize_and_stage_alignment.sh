#!/bin/bash

# USAGE:
# ./organize_and_stage_alignment.sh MOUSE_ID GSM_ID
# Example:
# ./organize_and_stage_alignment.sh mouse1 GSM8697658

MOUSE_ID="$1"
GSM_ID="$2"

SRC_DIR=~/BCR_IGE/BCR_IGE/MICE_sequences/${GSM_ID}_${MOUSE_ID}
DST_DIR=~/BCR_IGE/BCR_IGE/MICE_AA_Alignments

# Ensure destination exists
mkdir -p "$DST_DIR"

# Define filenames
AA_FASTA="${SRC_DIR}/full_variable_region_aa_${MOUSE_ID}.fasta"
RENAMED_AA_FASTA="${DST_DIR}/AA_${GSM_ID}.fasta"

# Check if AA fasta exists
if [ ! -f "$AA_FASTA" ]; then
    echo "❌ AA FASTA file not found at: $AA_FASTA"
    exit 1
fi

# Move and rename
mv "$AA_FASTA" "$RENAMED_AA_FASTA"

echo "✅ Renamed and moved:"
echo "   $AA_FASTA → $RENAMED_AA_FASTA"

# Optional: Suggest MAFFT alignment for local processing
echo ""
echo "📢 To align this with your master file on your Mac:"
echo "   cat all_sequences_combined_aligned.fasta AA_${GSM_ID}.fasta > tmp_combined.fasta"
echo "   mafft tmp_combined.fasta > all_sequences_combined_aligned_UPDATED.fasta"
echo "   rm tmp_combined.fasta"

echo ""
echo "📁 Done. Your file is now ready in:"
echo "   $DST_DIR"
