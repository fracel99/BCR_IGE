This project aims to find conserved and variable regions in the VDJ genes of mouse and/or human by analizing available data on BCRs, more specifically IgEs.
This GitHub repository is meant for storage, curation and organization of that data.
# Descriptions:
Commands: these are shell and python scripts to extract VDJ sequences from bigger files (e.g. CSV) and convert them to amino acid sequences.

MICE_Sequeneces: directory for curated DNA sequences.

MICE_AA_Sequences: the DNA sequences are tranlsated to AA sequences using biopython and then aligned with MAFFT. The result is uploaded to this folder and can be visualized with other programs such as Jalview.

Visualization: these files can be opened using Jalview and ChimeraX to see the conserved and variable regions across the MSA created from the sequences in "MICE_AA_Sequences".
