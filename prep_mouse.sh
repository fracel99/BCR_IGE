#!/bin/bash

# USAGE:
# ./prep_mouse.sh GEO_ID mouseX DOWNLOAD_LINK

# Read input arguments
GEO_ID="$1"
MOUSE_ID="$2"
DOWNLOAD_LINK="$3"
WORKDIR="${GEO_ID}_${MOUSE_ID}"

# Create working directory
cd ~/BCR/IGE/MICE || exit 1
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

# Download and extract the CSV file
echo "🔽 Downloading annotations..."
wget "$DOWNLOAD_LINK" -O "${WORKDIR}_filtered_contig_annotations.csv.gz"
gunzip "${WORKDIR}_filtered_contig_annotations.csv.gz"

# Filter only IGHE rows and store in separate file
echo "🧪 Filtering for IGHE..."
INPUT_CSV="${WORKDIR}_filtered_contig_annotations.csv"
IGE_CSV="${WORKDIR}_IgE.csv"
(head -n 1 "$INPUT_CSV" && grep IGHE "$INPUT_CSV") > "$IGE_CSV"

# Store IgE sequence count and V, D, J gene usage summary
echo "📊 Writing summary to ${MOUSE_ID}_Info.txt..."
echo "Number of IGHE sequences:" > "${MOUSE_ID}_Info.txt"
grep -c IGHE "$IGE_CSV" >> "${MOUSE_ID}_Info.txt"

echo -e "\nTop V genes used:" >> "${MOUSE_ID}_Info.txt"
cut -d',' -f7 "$IGE_CSV" | tail -n +2 | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

echo -e "\nTop D genes used:" >> "${MOUSE_ID}_Info.txt"
cut -d',' -f8 "$IGE_CSV" | tail -n +2 | grep -v '^$' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

echo -e "\nTop J genes used:" >> "${MOUSE_ID}_Info.txt"
cut -d',' -f9 "$IGE_CSV" | tail -n +2 | grep -v '^$' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

# Extract V region NTs (FWR1–FWR4 + CDR1–CDR3)
echo "🧬 Extracting full V region (NT)..."
cut -d',' -f14,16,18,20,22,24,26 "$IGE_CSV" | tail -n +2 > full_variable_region_nt_${MOUSE_ID}.txt
sed 's/,//g' full_variable_region_nt_${MOUSE_ID}.txt > full_variable_region_nt_${MOUSE_ID}_clean.txt

# Convert to FASTA format
awk '{print ">seq" NR "\n" $0}' full_variable_region_nt_${MOUSE_ID}_clean.txt > full_variable_region_nt_${MOUSE_ID}.fasta

# Clean up intermediate files
rm full_variable_region_nt_${MOUSE_ID}.txt
rm full_variable_region_nt_${MOUSE_ID}_clean.txt
rm "$IGE_CSV"

# Final note
echo "✅ Done: FASTA saved as full_variable_region_nt_${MOUSE_ID}.fasta"
echo "📁 Output directory: ~/BCR/IGE/MICE/${WORKDIR}"
