from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_file = "cdr3_nt_unknown.fasta"
output_file = "cdr3_aa_unknown.fasta"

translated_records = []

for record in SeqIO.parse(input_file, "fasta"):
    nt_seq = str(record.seq).replace("-", "")
    aa_seq = Seq(nt_seq).translate(to_stop=False)
    translated_record = SeqRecord(
        aa_seq,
        id=record.id,
        description=""
    )
    translated_records.append(translated_record)

SeqIO.write(translated_records, output_file, "fasta")
print(f"✅ Translated CDR3 sequences saved to {output_file}")
