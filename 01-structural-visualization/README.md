# Trastuzumab–HER2 Binding Interface Visualization

Structural visualization of the therapeutic antibody **trastuzumab** (Herceptin) bound to the extracellular domain of its target, the **HER2 receptor**, using PyMOL.

## Background
Trastuzumab is a monoclonal antibody used to treat HER2-positive breast cancer. This visualization highlights the CDR-H3 loop of the antibody heavy chain — the variable region most directly responsible for antigen recognition and a known hotspot for biopharmaceutical liabilities such as aggregation propensity and chemical degradation.

## Source
- PDB ID: 1N8Z
- Reference: Cho et al., *Nature* 2003

## Method
PyMOL session using open-source PyMOL 3.1.0. Chains assigned as:
- Chain A: trastuzumab heavy chain (Fab)
- Chain B: trastuzumab light chain (Fab)
- Chain C: HER2 receptor extracellular domain

CDR-H3 loop highlighted at residues 95–102 (heavy chain).

## Reproduce
\`\`\`bash
pymol visualize_trastuzumab.pml
\`\`\`

## Files
- `visualize_trastuzumab.pml` — PyMOL script
- `trastuzumab_HER2_interface.png` — rendered image
- `trastuzumab_session.pse` — saved PyMOL session
