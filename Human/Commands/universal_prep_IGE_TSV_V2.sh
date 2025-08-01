from pathlib import Path
from shutil import copyfile

# Redefine paths after kernel reset
original_script = Path("/mnt/data/universal_prep_IGE_HUMAN_V10.sh")
tsv_script = Path("/mnt/data/universal_prep_IGE_HUMAN_TSV.sh")

# Copy the original CSV-based script
copyfile(original_script, tsv_script)

# Read the content and replace all necessary parts
content = tsv_script.read_text()

# Replace .csv with .tsv where relevant
content = content.replace(".csv.gz", ".tsv.gz") \
                 .replace(".csv", ".tsv") \
                 .replace("filtered_contig_annotations.csv", "filtered_contig_annotations.tsv")

# Remove the "TSV to CSV" block entirely since the data is already in TSV
import re
content = re.sub(
    r'echo "🔁 Converting TSV to CSV \(if needed\)\.\.\."\n.*?fi\n\n', 
    '', content, flags=re.DOTALL
)

# Replace all comma delimiters (e.g., -F',' and cut -d',') with tab versions
content = content.replace("-F','", "-F'\t'") \
                 .replace("cut -d','", "cut -f") \
                 .replace("awk -F','", "awk -F'\t'")

# Save the modified version
tsv_script.write_text(content)

tsv_script.name  # Return the new script filename for confirmation
