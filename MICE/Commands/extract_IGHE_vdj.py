#!/usr/bin/env python3

import csv
import sys

# USAGE:
# ./extract_IGHE_vdj.py INPUT.tsv.gz OUTPUT.fasta MOUSE_ID

import gzip

input_file = sys.argv[1]
output_file = sys.argv[2]
mouse_id = sys.argv[3]

with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    count = 1
    for row in reader:
        c_call = row['c_call']
        if c_call != 'IGHE':
            continue
        sequence = row['sequence']
        try:
            v_start = int(row['v_sequence_start'])
            j_end = int(row['j_sequence_end'])
            vdj_seq = sequence[v_start:j_end]
        except (ValueError, IndexError):
            continue  # skip problematic rows
        if len(vdj_seq) < 100:
            continue  # too short to be real
        header = f">{mouse_id}_seq{count}"
        outfile.write(f"{header}\n{vdj_seq}\n")
        count += 1
