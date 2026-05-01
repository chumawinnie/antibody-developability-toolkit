# PyMOL script: Visualize trastuzumab Fab bound to HER2 receptor
# PDB: 1N8Z (Cho et al., 2003)
# Highlights CDR-H3 loop at the antibody-antigen interface

fetch 1n8z
remove resn HOH

# Chain A = heavy chain, B = light chain, C = HER2 ECD
color salmon, chain A
color skyblue, chain B
color palegreen, chain C
bg_color white

# Highlight CDR-H3 binding loop
color red, chain A and resi 95-102
show sticks, chain A and resi 95-102

# Frame the binding interface
orient chain A and resi 95-102
zoom chain A and resi 95-102, 15

# Render publication-quality image
ray 1600, 1200
png trastuzumab_HER2_interface.png
