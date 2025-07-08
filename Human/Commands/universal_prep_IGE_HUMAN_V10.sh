#!/bin/bash

# USAGE:
# ./universal_prep_IGE_HUMAN_V5.sh GEO_ID HUMAN_ID DOWNLOAD_LINK

GEO_ID="$1"
HUMAN_ID="$2"
DOWNLOAD_LINK="$3"
WORKDIR="${GEO_ID}_${HUMAN_ID}"

cd ~/BCR_IGE/BCR_IGE/Human/DNA_sequences || exit 1
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "🔽 Downloading file..."
wget "$DOWNLOAD_LINK" -O "${WORKDIR}_raw.csv.gz"

if ! gunzip "${WORKDIR}_raw.csv.gz"; then
    echo "❌ Failed to unzip file. Make sure it's a valid .csv.gz file. Exiting."
    exit 1
fi

CSV="${WORKDIR}_raw.csv"
mv "$CSV" "${WORKDIR}_filtered_contig_annotations.csv"
CSV="${WORKDIR}_filtered_contig_annotations.csv"

echo "🔁 Converting TSV to CSV (if needed)..."
DELIM=$(head -n 1 "$CSV" | grep -o $'\t' | wc -l)
if [ "$DELIM" -gt 0 ]; then
    sed 's/\t/,/g' "$CSV" > "${CSV%.csv}_converted.csv"
    CSV="${CSV%.csv}_converted.csv"
fi

echo "🧪 Filtering for IGHE..."
(head -n 1 "$CSV" && awk -F',' '$10 ~ /IGHE/' "$CSV") > "${WORKDIR}_IgE.csv"
IGE_CSV="${WORKDIR}_IgE.csv"

HEADER=$(head -n 1 "$IGE_CSV")
if echo "$HEADER" | grep -q "cdr3_nt"; then
    FORMAT="10x"
elif echo "$HEADER" | grep -q "IR_VDJ_1_junction"; then
    FORMAT="AIRR"
elif echo "$HEADER" | grep -q "v_call" && echo "$HEADER" | grep -q "cdr3"; then
    FORMAT="IGBLAST"
else
    echo "❌ Unknown CSV format. Exiting."
    exit 1
fi

echo "📊 Writing ${HUMAN_ID}_Info.txt..."
echo "Number of IGHE sequences:" > "${HUMAN_ID}_Info.txt"
grep -c IGHE "$IGE_CSV" >> "${HUMAN_ID}_Info.txt"

if [ "$FORMAT" = "10x" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f7 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f8 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f9 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo "📦 Generating IGHE_variable_region_nt_${HUMAN_ID}.fasta with barcodes..."
    tail -n +2 "$IGE_CSV" | awk -F',' '{
        seq = ""; 
        for (i=14; i<=26; i+=2) { 
            gsub(/"/, "", $i); seq = seq $i 
        }
        if (seq != "") print ">" $1 " | IGHE\n" seq
    }' > IGHE_variable_region_nt_${HUMAN_ID}.fasta

elif [ "$FORMAT" = "AIRR" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f28 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f32 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f36 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    cut -d',' -f16 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v 'None' > temp_cdr3.txt
    awk '{print ">seq" NR "\n" $0}' temp_cdr3.txt > cdr3_nt_${HUMAN_ID}.fasta
    rm temp_cdr3.txt

elif [ "$FORMAT" = "IGBLAST" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f5 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f6 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -d',' -f7 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    cut -d',' -f27 "$IGE_CSV" | tail -n +2 | grep -v '^$' > temp_seq.txt
    awk '{print ">seq" NR "\n" $0}' temp_seq.txt > IGHE_variable_region_nt_${HUMAN_ID}.fasta
    rm temp_seq.txt
fi

rm "$IGE_CSV"
echo "✅ Finished. Output in ~/BCR_IGE/BCR_IGE/Human/DNA_sequences/${WORKDIR}"

# --- IGK/IGL BLOCK ---
CSV="${WORKDIR}_filtered_contig_annotations.csv"
OUTPUT="IGK_IGL_variable_region_nt_${HUMAN_ID}.fasta"

echo "🔍 Extracting IGHE barcodes..."
awk -F',' 'NR>1 && $10 ~ /IGHE/ {print $1}' "$CSV" | sort | uniq > ige_barcodes.txt

echo "🔗 Matching IGK/IGL sequences with IGHE barcodes..."
awk -F',' 'NR==FNR {ige[$1]; next} 
    ($6 == "IGK" || $6 == "IGL") && ($1 in ige) {
        seq = ""; 
        for (i=14; i<=26; i+=2) if ($i != "") seq = seq $i;
        if (seq != "") print ">" $1 " | " $6 "\n" seq
    }' ige_barcodes.txt "$CSV" > "$OUTPUT"

echo "✅ IGK/IGL output: $OUTPUT"

# --- SPLIT IGK/IGL FILES ---
echo "✂️ Splitting IGK_IGL_variable_region_nt_${HUMAN_ID}.fasta..."

awk '/^>/{header=$0; getline seq; if (header ~ /IGK$/) {print header "\n" seq >> "IGK_variable_region_nt_'${HUMAN_ID}'.fasta"} else if (header ~ /IGL$/) {print header "\n" seq >> "IGL_variable_region_nt_'${HUMAN_ID}'.fasta"}}' \
    IGK_IGL_variable_region_nt_${HUMAN_ID}.fasta

echo "✅ IGK and IGL FASTA files created:"
echo " - IGK_variable_region_nt_${HUMAN_ID}.fasta"
echo " - IGL_variable_region_nt_${HUMAN_ID}.fasta"

# --- LIGHT CHAIN GENE USAGE SUMMARY (IGK/IGL matched to IGHE) ---
# 🧬 Separate IGK and IGL gene usage only for IGHE-paired barcodes
echo -e "\nTop V genes used by IGK paired with IGHE:" >> "${HUMAN_ID}_Info.txt"
awk -F',' 'NR>1 && $10 == "IGHE" {seen[$1]++} END {for (i in seen) print i}' "$CSV" > ige_barcodes.txt
awk -F',' 'NR==FNR {barcodes[$1]; next} $6 == "IGK" && ($1 in barcodes) {gsub(/"/, "", $7); print $7}' ige_barcodes.txt "$CSV" \
  | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

echo -e "\nTop J genes used by IGK paired with IGHE:" >> "${HUMAN_ID}_Info.txt"
awk -F',' 'NR==FNR {barcodes[$1]; next} $6 == "IGK" && ($1 in barcodes) {gsub(/"/, "", $9); print $9}' ige_barcodes.txt "$CSV" \
  | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

echo -e "\nTop V genes used by IGL paired with IGHE:" >> "${HUMAN_ID}_Info.txt"
awk -F',' 'NR==FNR {barcodes[$1]; next} $6 == "IGL" && ($1 in barcodes) {gsub(/"/, "", $7); print $7}' ige_barcodes.txt "$CSV" \
  | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

echo -e "\nTop J genes used by IGL paired with IGHE:" >> "${HUMAN_ID}_Info.txt"
awk -F',' 'NR==FNR {barcodes[$1]; next} $6 == "IGL" && ($1 in barcodes) {gsub(/"/, "", $9); print $9}' ige_barcodes.txt "$CSV" \
  | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

# --- FINAL PAIRING BLOCK ---
echo "🔗 Pairing IGHE and IGK/IGL sequences by barcode..."

IGH_FILE="IGHE_variable_region_nt_${HUMAN_ID}.fasta"
IGKL_FILE="IGK_IGL_variable_region_nt_${HUMAN_ID}.fasta"
OUTPUT_PAIRS="IGE_variable_region_pairs_nt_${HUMAN_ID}.fasta"

# Convert to TSV: barcode \t sequence
awk '/^>/{sub(/^>/,""); split($0,a,/[\t ]/); header=a[1]; next} {print header "\t" $0}' "$IGH_FILE" > IGHE_seqs.tsv
awk '/^>/{sub(/^>/,""); split($0,a,/[\t| ]/); header=a[1]; type=$NF; next} {print header "\t" $0 "\t" type}' "$IGKL_FILE" > IGKL_seqs.tsv

# Sort both files by barcode
sort IGHE_seqs.tsv > IGHE_sorted.tsv
sort IGKL_seqs.tsv > IGKL_sorted.tsv

# Join and format as FASTA
join -t $'\t' IGHE_sorted.tsv IGKL_sorted.tsv | awk -F'\t' '{
    print ">" $1 " IGH\n" $2 "\n>" $1 " " $4 "\n" $3
}' > "$OUTPUT_PAIRS"

# Cleanup
rm IGHE_seqs.tsv IGKL_seqs.tsv IGHE_sorted.tsv IGKL_sorted.tsv

echo "✅ Paired output written to: $OUTPUT_PAIRS"

# --- TRANSLATE TO AMINO ACIDS ---
echo -e "\n🧬 Translating nucleotide sequences to amino acids..."

TRANSLATE_SCRIPT="${HOME}/BCR_IGE/BCR_IGE/Human/Commands/translate_fasta.py"
TRANSLATE_DIR="${PWD}"  # current dataset directory

for CHAIN in IGHE IGK IGL; do
    INPUT="${TRANSLATE_DIR}/${CHAIN}_variable_region_nt_${HUMAN_ID}.fasta"
    OUTPUT="${TRANSLATE_DIR}/${CHAIN}_variable_region_aa_${HUMAN_ID}.fasta"

    if [ -f "$INPUT" ]; then
        python3 "$TRANSLATE_SCRIPT" "$INPUT" "$OUTPUT"
    else
        echo "⚠️ Skipping $CHAIN – no input file found: $INPUT"
    fi
done

echo "✅ Translation complete:"
echo "  - IGHE_variable_region_aa_${HUMAN_ID}.fasta"
echo "  - IGK_variable_region_aa_${HUMAN_ID}.fasta"
echo "  - IGL_variable_region_aa_${HUMAN_ID}.fasta"
