#!/bin/bash

# USAGE:
# ./universal_prep_IGE_TSV.sh GEO_ID HUMAN_ID DOWNLOAD_LINK

GEO_ID="$1"
HUMAN_ID="$2"
DOWNLOAD_LINK="$3"
WORKDIR="${GEO_ID}_${HUMAN_ID}"

cd ~/BCR_IGE/BCR_IGE/Human/DNA_sequences || exit 1
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "🔽 Downloading file..."
FILENAME="${WORKDIR}_raw.tsv.gz"
wget "$DOWNLOAD_LINK" -O "$FILENAME"

if ! gunzip "$FILENAME"; then
    echo "❌ Failed to unzip file. Exiting."
    exit 1
fi

TSV="${FILENAME%.gz}"
mv "$TSV" "${WORKDIR}_filtered_contig_annotations.tsv"
TSV="${WORKDIR}_filtered_contig_annotations.tsv"

echo "🧪 Filtering for IGHE..."
(head -n 1 "$TSV" && awk -F'\t' '$14 ~ /IGHE/' "$TSV") > "${WORKDIR}_IgE.tsv"
IGE_TSV="${WORKDIR}_IgE.tsv"

HEADER=$(head -n 1 "$IGE_TSV")
if echo "$HEADER" | grep -q "cdr3_nt"; then
    FORMAT="10x"
elif echo "$HEADER" | grep -q "junction_aa" && echo "$HEADER" | grep -q "v_call"; then
    FORMAT="AIRR"
elif echo "$HEADER" | grep -q "v_call" && echo "$HEADER" | grep -q "cdr3"; then
    FORMAT="IGBLAST"
else
    echo "❌ Unknown TSV format. Exiting."
    exit 1
fi

echo "📊 Writing ${HUMAN_ID}_Info.txt..."
echo "Number of IGHE sequences:" > "${HUMAN_ID}_Info.txt"
grep -c IGHE "$IGE_TSV" >> "${HUMAN_ID}_Info.txt"

if [ "$FORMAT" = "AIRR" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f28 "$IGE_TSV" | tail -n +2 | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f32 "$IGE_TSV" | tail -n +2 | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f36 "$IGE_TSV" | tail -n +2 | grep -v '^$' | grep -v nan | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    cut -f16 "$IGE_TSV" | tail -n +2 | grep -v '^$' | grep -v 'None' > temp_cdr3.txt
    awk '{print ">seq" NR "\n" $0}' temp_cdr3.txt > cdr3_nt_${HUMAN_ID}.fasta
    rm temp_cdr3.txt
else
    echo "⚠️ Current TSV script only supports AIRR format."
fi

rm "$IGE_TSV"
echo "✅ Finished. Output in ~/BCR_IGE/BCR_IGE/Human/DNA_sequences/${WORKDIR}"
