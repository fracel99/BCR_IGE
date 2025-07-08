#!/bin/bash

# USAGE:
# ./universal_prep_IGE_TSV_V3.sh GEO_ID HUMAN_ID DOWNLOAD_LINK

GEO_ID="$1"
HUMAN_ID="$2"
DOWNLOAD_LINK="$3"
WORKDIR="${GEO_ID}_${HUMAN_ID}"

cd ~/BCR_IGE/BCR_IGE/Human/DNA_sequences || exit 1
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "🔽 Downloading file..."
wget "$DOWNLOAD_LINK" -O "${WORKDIR}_raw.tsv.gz"

if ! gunzip "${WORKDIR}_raw.tsv.gz"; then
    echo "❌ Failed to unzip file. Exiting."
    exit 1
fi

TSV="${WORKDIR}_raw.tsv"
mv "$TSV" "${WORKDIR}_filtered_contig_annotations.tsv"
TSV="${WORKDIR}_filtered_contig_annotations.tsv"

echo "🧪 Filtering for IGHE..."
(head -n 1 "$TSV" && awk -F'\t' '$14 ~ /IGHE/' "$TSV") > "${WORKDIR}_IgE.csv"
IGE_TSV="${WORKDIR}_IgE.tsv"

HEADER=$(head -n 1 "$IGE_TSV")
if echo "$HEADER" | grep -q $'\tcdr3_nt'; then
    FORMAT="10x"
elif echo "$HEADER" | grep -q $'\tIR_VDJ_1_junction'; then
    FORMAT="AIRR"
elif echo "$HEADER" | grep -q $'\tv_call' && echo "$HEADER" | grep -q $'\tjunction'; then
    FORMAT="AIRR_TSV"
elif echo "$HEADER" | grep -q $'\tv_call' && echo "$HEADER" | grep -q $'\tcdr3'; then
    FORMAT="IGBLAST"
else
    echo "❌ Unknown TSV format. Exiting."
    exit 1
fi

echo "📊 Writing ${HUMAN_ID}_Info.txt..."
echo "Number of IGHE sequences:" > "${HUMAN_ID}_Info.txt"
grep -c IGHE "$IGE_TSV" >> "${HUMAN_ID}_Info.txt"

if [ "$FORMAT" = "10x" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f7 "$IGE_TSV" | tail -n +2 | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f8 "$IGE_TSV" | tail -n +2 | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f9 "$IGE_TSV" | tail -n +2 | grep -v '^$' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo "📦 Generating IGHE_variable_region_nt_${HUMAN_ID}.fasta..."
    tail -n +2 "$IGE_TSV" | awk -F'\t' '{
        seq = ""; 
        for (i=14; i<=26; i+=2) { seq = seq $i }
        if (seq != "") print ">" $1 " | IGHE\n" seq
    }' > IGHE_variable_region_nt_${HUMAN_ID}.fasta

elif [ "$FORMAT" = "AIRR" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f28 "$IGE_TSV" | tail -n +2 | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f32 "$IGE_TSV" | tail -n +2 | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f36 "$IGE_TSV" | tail -n +2 | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    cut -f16 "$IGE_TSV" | tail -n +2 | grep -v '^$' > temp_cdr3.txt
    awk '{print ">seq" NR "\n" $0}' temp_cdr3.txt > cdr3_nt_${HUMAN_ID}.fasta
    rm temp_cdr3.txt

elif [ "$FORMAT" = "IGBLAST" ]; then
    echo -e "\nTop V genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f5 "$IGE_TSV" | tail -n +2 | sed 's/ .*//' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop D genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f6 "$IGE_TSV" | tail -n +2 | sed 's/ .*//' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    echo -e "\nTop J genes used:" >> "${HUMAN_ID}_Info.txt"
    cut -f7 "$IGE_TSV" | tail -n +2 | sed 's/ .*//' | sort | uniq -c | sort -nr >> "${HUMAN_ID}_Info.txt"

    cut -f27 "$IGE_TSV" | tail -n +2 | grep -v '^$' > temp_seq.txt
    awk '{print ">seq" NR "\n" $0}' temp_seq.txt > IGHE_variable_region_nt_${HUMAN_ID}.fasta
    rm temp_seq.txt
fi

rm "$IGE_TSV"
echo "✅ Finished. Output in ~/BCR_IGE/BCR_IGE/Human/DNA_sequences/${WORKDIR}"
