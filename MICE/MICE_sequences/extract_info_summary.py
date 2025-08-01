import pandas as pd
import sys
from collections import Counter

input_file = sys.argv[1]
df = pd.read_csv(input_file, sep="\t")

# Drop rows with missing c_call
df = df.dropna(subset=["c_call"])

# Filter IGHE
df = df[df["c_call"].str.upper().str.startswith("IGHE")]

# 🔽 STRIP ALLELE INFO (*01, *02, etc.) — insert here
df["v_call"] = df["v_call"].str.replace(r"\*.*", "", regex=True)
df["d_call"] = df["d_call"].str.replace(r"\*.*", "", regex=True)
df["j_call"] = df["j_call"].str.replace(r"\*.*", "", regex=True)

# Count total IGHE
total_ighe = len(df)

# Count genes (handle comma-separated multi-assignments safely)
def flatten_column(col):
    return [item.strip() for entry in col.dropna() for item in str(entry).split(",")]

v_counts = Counter(flatten_column(df["v_call"]))
d_counts = Counter(flatten_column(df["d_call"]))
j_counts = Counter(flatten_column(df["j_call"]))

# Create summary output
with open("VDJ_summary.txt", "w") as out:
    out.write(f"Number of IGHE sequences: {total_ighe}\n\n")

    out.write("Top V genes used:\n")
    for gene, count in v_counts.most_common():
        out.write(f"    {count:>3} {gene}\n")

    out.write("\nTop D genes used:\n")
    for gene, count in d_counts.most_common():
        out.write(f"    {count:>3} {gene}\n")

    out.write("\nTop J genes used:\n")
    for gene, count in j_counts.most_common():
        out.write(f"    {count:>3} {gene}\n")
