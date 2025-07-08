import pandas as pd
import sys
from pathlib import Path

def process_file(input_path, mouse_id):
    df = pd.read_csv(input_path)

    # Flexible C gene column detection
    c_call_column = None
    for option in ['C_CALL', 'C_CALL.H', 'C_CALL.L', 'c_gene']:
        if option in df.columns:
            c_call_column = option
            break

    if not c_call_column:
        print("No C_CALL column found.")
        return

    # Filter IGHE
    ige_df = df[df[c_call_column].str.contains('IGHE', na=False, case=False)]

    if ige_df.empty:
        print("No IGHE sequences found.")
        return

    output_dir = Path(input_path).parent

    # Write FASTA
    fasta_path = output_dir / f"full_variable_region_nt_{mouse_id}.fasta"
    with open(fasta_path, 'w') as f:
        for i, row in ige_df.iterrows():
            sequence = row.get('SEQUENCE_INPUT', row.get('cdr3_nt', ''))
            f.write(f">{mouse_id}_seq{i}\n{sequence}\n")

    # Write Info
    info_path = output_dir / f"{mouse_id}_Info.txt"
    with open(info_path, 'w') as f:
        f.write(f"Number of IGHE sequences:\n{len(ige_df)}\n\n")
        for gene in ['v_gene', 'd_gene', 'j_gene']:
            if gene in df.columns:
                f.write(f"Top {gene} usage:\n")
                f.write(ige_df[gene].value_counts().to_string())
                f.write("\n\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_local_ige_v4.py <input_csv> <mouse_id>")
    else:
        process_file(sys.argv[1], sys.argv[2])
