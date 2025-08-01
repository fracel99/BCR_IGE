#!/bin/bash

# USAGE:
# ./download_TAR.sh <TAR_DOWNLOAD_LINK>

# Example:
# ./download_TAR.sh https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193868/suppl/GSE193868_A_SPL.tar.gz

TAR_URL="$1"
if [[ -z "$TAR_URL" ]]; then
    echo "❌ Please provide a TAR download link."
    exit 1
fi

FILENAME=$(basename "$TAR_URL")                   # e.g. GSE193868_A_SPL.tar.gz
GEO_ID=$(echo "$FILENAME" | cut -d'_' -f1)         # e.g. GSE193868
SUBID_PART=$(echo "$FILENAME" | cut -d'_' -f2)

# Decide SUBID
if [[ "$SUBID_PART" == *.tar.gz || -z "$SUBID_PART" ]]; then
    SUBID="$GEO_ID"
else
    SUBID=$(echo "$FILENAME" | sed -E "s/^${GEO_ID}_//" | sed 's/.tar.gz$//')
fi

WORKDIR="${HOME}/BCR_IGE/BCR_IGE/Human/DNA_sequences/${GEO_ID}/${SUBID}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "�� Downloading $FILENAME ..."
wget "$TAR_URL" -O "$FILENAME"

echo "📦 Extracting $FILENAME ..."
tar -xzf "$FILENAME"

echo "🔍 Searching for .csv annotation file..."
ANNOTATION_FILE=$(find . -type f -name "*filtered_contig*.csv" | head -n 1)

if [[ -z "$ANNOTATION_FILE" ]]; then
    echo "❌ No CSV file found in archive. Exiting."
    exit 1
fi

echo "✅ Found annotation file: $ANNOTATION_FILE"
echo "🚀 Running universal_prep_IGE_HUMAN_V11.sh ..."
~/BCR_IGE/BCR_IGE/Human/Commands/universal_prep_IGE_HUMAN_V11.sh "$GEO_ID" "$SUBID" "$(realpath "$ANNOTATION_FILE")"
