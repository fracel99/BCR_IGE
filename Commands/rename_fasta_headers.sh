#!/bin/bash

# USAGE:
# ./rename_fasta_headers.sh INPUT.fasta MOUSE_ID OUTPUT.fasta

INPUT="$1"
MOUSE_ID="$2"
OUTPUT="$3"

awk -v id="$MOUSE_ID" '
  /^>/ {
    i++
    print ">" id "_seq" i
    next
  }
  { print }
' "$INPUT" > "$OUTPUT"
