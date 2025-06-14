#!/bin/bash

# USAGE:
# ./universal_prep_IGE_V4.sh GEO_ID mouseX DOWNLOAD_LINK

GEO_ID="$1"
MOUSE_ID="$2"
DOWNLOAD_LINK="$3"
WORKDIR="${GEO_ID}_${MOUSE_ID}"

cd ~/BCR_IGE/BCR_IGE/MICE_sequences || exit 1
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
(head -n 1 "$CSV" && grep IGHE "$CSV") > "${WORKDIR}_IgE.csv"
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

echo "📊 Writing ${MOUSE_ID}_Info.txt..."
echo "Number of IGHE sequences:" > "${MOUSE_ID}_Info.txt"
grep -c IGHE "$IGE_CSV" >> "${MOUSE_ID}_Info.txt"

if [ "$FORMAT" = "10x" ]; then
    echo -e "\nTop V genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f7 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f8 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f9 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    cut -d',' -f14,16,18,20,22,24,26 "$IGE_CSV" | tail -n +2 > temp_nt.txt
    sed 's/,//g' temp_nt.txt > temp_clean.txt
    awk '{print ">seq" NR "\n" $0}' temp_clean.txt > full_variable_region_nt_${MOUSE_ID}.fasta
    rm temp_nt.txt temp_clean.txt

elif [ "$FORMAT" = "AIRR" ]; then
    echo -e "\nTop V genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f28 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f32 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f36 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v nan | sed 's/"//g' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    cut -d',' -f16 "$IGE_CSV" | tail -n +2 | grep -v '^$' | grep -v 'None' > temp_cdr3.txt
    awk '{print ">seq" NR "\n" $0}' temp_cdr3.txt > cdr3_nt_${MOUSE_ID}.fasta
    rm temp_cdr3.txt

elif [ "$FORMAT" = "IGBLAST" ]; then
    echo -e "\nTop V genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f5 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f6 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${MOUSE_ID}_Info.txt"
    cut -d',' -f7 "$IGE_CSV" | tail -n +2 | sed 's/ .*//' | grep -v '^$' | sort | uniq -c | sort -nr >> "${MOUSE_ID}_Info.txt"

    cut -d',' -f27 "$IGE_CSV" | tail -n +2 | grep -v '^$' > temp_seq.txt
    awk '{print ">seq" NR "\n" $0}' temp_seq.txt > full_variable_region_nt_${MOUSE_ID}.fasta
    rm temp_seq.txt
fi

rm "$IGE_CSV"
echo "✅ Finished. Output in ~/BCR_IGE/BCR_IGE/MICE_sequences/${WORKDIR}"
