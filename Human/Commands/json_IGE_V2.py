import json
import sys
from collections import defaultdict, Counter

# --- Load JSON ---
with open(sys.argv[1]) as f:
    data = json.load(f)

# --- Organize annotations by barcode ---
barcodes = defaultdict(list)

for record in data:
    if not record.get("is_cell") or not record.get("productive"):
        continue
    barcode = record["barcode"]
    annotations = record.get("annotations", [])
    chains = defaultdict(list)

    for ann in annotations:
        feat = ann.get("feature", {})
        chain = feat.get("chain")
        region = feat.get("region_type")
        gene = feat.get("gene_name")

        if chain and region and gene:
            if region.endswith("V-REGION"):
                chains[chain].append(("V", gene))
            elif region == "D-REGION":
                chains[chain].append(("D", gene))
            elif region == "J-REGION":
                chains[chain].append(("J", gene))

    for chain, regions in chains.items():
        barcodes[barcode].append((chain, dict(regions)))

# --- Pair heavy and light chains ---
pairings = []
gene_counts = {
    "IGH_V": Counter(),
    "IGH_D": Counter(),
    "IGH_J": Counter(),
    "IGK_V": Counter(),
    "IGK_J": Counter(),
    "IGL_V": Counter(),
    "IGL_J": Counter(),
}
ighe_count = 0

for barcode, chains in barcodes.items():
    chain_dict = {c[0]: c[1] for c in chains}

    if "IGH" in chain_dict and any(x in chain_dict for x in ["IGK", "IGL"]):
        hc = chain_dict["IGH"]
        c_gene = next((ann["feature"]["gene_name"] for rec in data if rec["barcode"] == barcode
                       for ann in rec["annotations"] if ann["feature"]["chain"] == "IGH" and
                       ann["feature"]["region_type"] == "C-REGION"), "")
        if c_gene.startswith("IGHE"):
            ighe_count += 1

            lc_type = "IGK" if "IGK" in chain_dict else "IGL"
            lc = chain_dict.get(lc_type, {})

            pairings.append([
                barcode,
                hc.get("V", ""),
                hc.get("D", ""),
                hc.get("J", ""),
                lc_type,
                lc.get("V", ""),
                lc.get("D", ""),
                lc.get("J", "")
            ])

            gene_counts["IGH_V"][hc.get("V", "")] += 1
            gene_counts["IGH_D"][hc.get("D", "")] += 1
            gene_counts["IGH_J"][hc.get("J", "")] += 1
            gene_counts[f"{lc_type}_V"][lc.get("V", "")] += 1
            gene_counts[f"{lc_type}_J"][lc.get("J", "")] += 1

# --- Write output files ---
with open("VDJ_pairings_ID.tsv", "w") as f:
    f.write("barcode\tIGH_V\tIGH_D\tIGH_J\tLight_chain\tLight_V\tLight_D\tLight_J\n")
    for row in pairings:
        f.write("\t".join(row) + "\n")

with open("VDJ_pairings_summary.tsv", "w") as f:
    counter = Counter([" | ".join(row[1:]) for row in pairings])
    for pair, count in counter.most_common():
        f.write(f"{count}\t{pair}\n")

with open("ID_info.txt", "w") as f:
    f.write(f"Number of IGHE sequences: {ighe_count}\n\n")
    f.write("Top V genes used:\n")
    for v, c in gene_counts["IGH_V"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop D genes used:\n")
    for v, c in gene_counts["IGH_D"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop J genes used:\n")
    for v, c in gene_counts["IGH_J"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop V genes used by IGK paired with IGHE:\n")
    for v, c in gene_counts["IGK_V"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop J genes used by IGK paired with IGHE:\n")
    for v, c in gene_counts["IGK_J"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop V genes used by IGL paired with IGHE:\n")
    for v, c in gene_counts["IGL_V"].most_common(5):
        f.write(f"{v}\t{c}\n")
    f.write("\nTop J genes used by IGL paired with IGHE:\n")
    for v, c in gene_counts["IGL_J"].most_common(5):
        f.write(f"{v}\t{c}\n")
