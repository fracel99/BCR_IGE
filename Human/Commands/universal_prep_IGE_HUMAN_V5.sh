#!/bin/bash

# USAGE:
# ./universal_prep_IGE_HUMAN.sh GEO_ID HUMAN_ID DOWNLOAD_LINK

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

# IGHE-specific FASTA (using barcode headers)
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
        if (seq != "") print ">" $1 "\n" seq
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
        if (seq != "") print ">" $1 "\n" seq
    }' ige_barcodes.txt "$CSV" > "$OUTPUT"

echo "✅ Done. Output: $OUTPUT"

# --- FINAL PAIRING BLOCK ---
echo "🔗 Pairing IGHE and IGK/IGL sequences by barcode..."
IGH_FILE="IGHE_variable_region_nt_${HUMAN_ID}.fasta"
IGKL_FILE="IGK_IGL_variable_region_nt_${HUMAN_ID}.fasta"
OUTPUT_PAIRS="IGE_variable_region_pairs_nt_${HUMAN_ID}.fasta"

awk '/^>/{header=$0; next} {print header "\t" $0}' "$IGH_FILE" > IGH_temp.tsv
awk '/^>/{header=$0; next} {print header "\t" $0}' "$IGKL_FILE" > IGKL_temp.tsv
sort IGH_temp.tsv > IGH_sorted.tsv
sort IGKL_temp.tsv > IGKL_sorted.tsv
join -t$'\t' -j 1 IGH_sorted.tsv IGKL_sorted.tsv > joined.tsv

awk -F'\t' '{
    gsub(/^>/, "", $1);
    print ">" $1 " IGH\n" $2 "\n>" $1 " IGK/L\n" $3
}' joined.tsv > "$OUTPUT_PAIRS"

rm IGH_temp.tsv IGKL_temp.tsv IGH_sorted.tsv IGKL_sorted.tsv joined.tsv

echo "✅ Paired output: $OUTPUT_PAIRS"

