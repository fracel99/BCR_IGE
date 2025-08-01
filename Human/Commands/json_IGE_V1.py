import json
import sys
from collections import Counter, defaultdict

# Input file
input_file = sys.argv[1]

# Load JSON
with open(input_file) as f:
    data = json.load(f)

# Prepare containers
ighe_entries = []
pairings = []

# Group contigs by barcode
barcode_dict = defaultdict(list)
for entry in data:
    if entry.get("high_confidence") and entry.get("productive") and entry.get("is_cell"):
        barcode_dict[entry["barcode"]].append(entry)

# Process barcode groups
for barcode, contigs in barcode_dict.items():
    ig_heavy = None
    ig_light = None

    for c in contigs:
        if c.get("chain") == "IGH" and c.get("c_gene", "").startswith("IGHE"):
            ig_heavy = c
        elif c.get["chain"] in {"IGK", "IGL"}:
            ig_light = c

    if ig_heavy:
        ighe_entries.append(ig_heavy)

        # Handle pairing if available
        light = ig_light if ig_light else {}
        pairings.append({
            "barcode": barcode,
            "IGH_V": ighe_entries[-1].get("v_gene", "NA"),
            "IGH_D": ighe_entries[-1].get("d_gene", "NA"),
            "IGH_J": ighe_entries[-1].get("j_gene", "NA"),
            "Light_chain": light.get("chain", "NA"),
            "Light_V": light.get("v_gene", "NA"),
            "Light_D": light.get("d_gene", "NA"),
            "Light_J": light.get("j_gene", "NA"),
        })

# Helper to count genes
def count_genes(entries, chain="IGH", gene="v_gene"):
    return Counter([e.get(gene) for e in entries if e.get("chain") == chain and e.get(gene)])

# ID_info.txt
with open("ID_info.txt", "w") as out:
    out.write(f"Number of IGHE sequences: {len(ighe_entries)}\n\n")

    out.write("Top V genes used:\n")
    for gene, count in count_genes(ighe_entries, "IGH", "v_gene").most_common():
        out.write(f"    {count} {gene}\n")

    out.write("\nTop D genes used:\n")
    for gene, count in count_genes(ighe_entries, "IGH", "d_gene").most_common():
        out.write(f"    {count} {gene}\n")

    out.write("\nTop J genes used:\n")
    for gene, count in count_genes(ighe_entries, "IGH", "j_gene").most_common():
        out.write(f"    {count} {gene}\n")

    # Light chains paired with IGHE
    for light_chain in ["IGK", "IGL"]:
        paired = [p for p in pairings if p["Light_chain"] == light_chain]
        out.write(f"\nTop V genes used by {light_chain} paired with IGHE:\n")
        for gene, count in Counter([p["Light_V"] for p in paired if p["Light_V"] != "NA"]).most_common():
            out.write(f"    {count} {gene}\n")
        out.write(f"\nTop J genes used by {light_chain} paired with IGHE:\n")
        for gene, count in Counter([p["Light_J"] for p in paired if p["Light_J"] != "NA"]).most_common():
            out.write(f"    {count} {gene}\n")

# VDJ_pairings_ID.tsv
with open("VDJ_pairings_ID.tsv", "w") as out:
    out.write("barcode\tIGH_V\tIGH_D\tIGH_J\tLight_chain\tLight_V\tLight_D\tLight_J\n")
    for p in pairings:
        out.write("\t".join([p[k] for k in ["barcode", "IGH_V", "IGH_D", "IGH_J", "Light_chain", "Light_V", "Light_D", "Light_J"]]) + "\n")

# VDJ_pairings_summary.tsv
pair_format = lambda p: f"{p['IGH_V']}+{p['IGH_D']}+{p['IGH_J']} | {p['Light_chain']}:{p['Light_V']}+{p['Light_D']}+{p['Light_J']}"
summary_counts = Counter([pair_format(p) for p in pairings])

with open("VDJ_pairings_summary.tsv", "w") as out:
    out.write("Count\tPair\n")
    for pair, count in summary_counts.most_common():
        out.write(f"{count}\t{pair}\n")
