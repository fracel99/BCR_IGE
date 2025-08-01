import pandas as pd
import sys
from pathlib import Path

def process_file(input_path, mouse_id):
    # Read the CSV file
    df = pd.read_csv(input_path)

    # Determine which C_CALL column exists
    if 'c_gene' in df.columns:
        c_call_column = 'c_gene'
    elif 'C_CALL' in df.columns:
        c_call_column = 'C_CALL'
    elif 'C_CALL.H' in df.columns:
        c_call_column = 'C_CALL.H'
    elif 'C_CALL.L' in df.columns:
        c_call_column = 'C_CALL.L'
    else:
        print("No C_CALL column found.")
        return

    # Filter for IgE entries
    ige_df = df[df[c_call_column].str.contains('IGHE', na=False, case=False)]

    if ige_df.empty:
        print("No IGHE sequences found.")
        return

    # Output directory
    output_dir = Path(input_path).parent

    # Write full variable region FASTA (skip NaNs)
    fasta_path = output_dir / f"full_variable_region_nt_{mouse_id}.fasta"
    with open(fasta_path, 'w') as f:
        for i, row in ige_df.iterrows():
            seq = row.get('cdr3_nt', '')
            if pd.notna(seq):
                f.write(f">{mouse_id}_seq{i}\n{seq}\n")

    # Write info file
    info_path = output_dir / f"{mouse_id}_Info.txt"
    with open(info_path, 'w') as f:
        f.write(f"Number of IGHE sequences:\n{len(ige_df)}\n\n")

        for gene in ['v_gene', 'd_gene', 'j_gene']:
            if gene in ige_df.columns:
                f.write(f"Top {gene} usage:\n")
                f.write(ige_df[gene].value_counts().to_string())
                f.write("\n\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_local_ige_v4.py <input_csv> <mouse_id>")
    else:
        process_file(sys.argv[1], sys.argv[2])
