# TP53-Structural-Analysis
Structural comparison of the p53 DNA-binding domain across wild-type and cancer-associated mutant structures, using R and bio3d.

What it does
Amino acid composition (WT vs R175H)
B-factor flexibility across structures, with L2/L3 loops highlighted
Delta B-factor: mutant vs WT
Per-residue RMSD (Kabsch superposition)
Pairwise RMSD heatmap
Hotspot mutation frequency chart


Output
Six PNG plots (composition, flexibility, delta B-factor, RMSD, RMSD heatmap, hotspot frequencies).

Requirements
rinstall.packages(c("bio3d", "tibble", "scales", "cowplot", "viridis", "ggplot2", "dplyr", "tidyr"))
