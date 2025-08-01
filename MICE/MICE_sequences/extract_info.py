import pandas as pd
import sys
from collections import Counter

# Input
input_file = sys.argv[1]

# Load tsv
df = pd.read_csv(input_file, sep="\t", dtype=str)

# Filter IGHE
df = df[df["c_call"].fillna('').str.upper().str.startswith("IGHE")]

# Drop rows with missing essential data
df = df.dropna(subset=["v_call", "d_call", "j_call", "junction", "sequence"])

# Count clones by identical VDJ + junction + sequence
group_cols = ["v_call", "d_call", "j_call", "junction", "sequence"]
df["clone_count"] = df.groupby(group_cols)["sequence"].transform("count")

# Drop duplicates (one row per unique clone)
df_unique = df.drop_duplicates(subset=group_cols)

# Write ID_Info.txt
df_unique[["clone_count", "v_call", "d_call", "j_call", "junction"]].to_csv("ID_Info.txt", sep="\t", index=False)

# Write IGHE_variable_region_nt_ID.fasta
with open("IGHE_variable_region_nt_ID.fasta", "w") as f_out:
    for _, row in df_unique.iterrows():
        header = f">{row['v_call']}+{row['d_call']}+{row['j_call']}_count{row['clone_count']}"
        sequence = row["sequence"]
        f_out.write(f"{header}\n{sequence}\n")
