import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv) != 3:
    print("Usage: python translate_fasta.py input_nt.fasta output_aa.fasta")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

translated_records = []

for record in SeqIO.parse(input_file, "fasta"):
    nt_seq = str(record.seq).replace("-", "")
    nt_seq = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)]  # trim to multiple of 3
    aa_seq = Seq(nt_seq).translate(to_stop=False)
    translated_record = SeqRecord(
        aa_seq,
        id=record.id,
        description=""
    )
    translated_records.append(translated_record)

SeqIO.write(translated_records, output_file, "fasta")
print(f"✅ Translated {len(translated_records)} sequences to {output_file}")
