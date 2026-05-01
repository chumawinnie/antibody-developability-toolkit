# Trastuzumab Biophysical Liability Analysis

Computational identification of formulation-relevant liabilities on the therapeutic antibody trastuzumab (PDB 1N8Z), using open-source Python tools.

## Method

Combines sequence chemistry with structural surface-exposure to flag three classes of biopharmaceutical liabilities:

1. Aggregation hotspots: surface-exposed hydrophobic residues (SASA at least 20 square Angstrom)
2. Deamidation sites: Asn-Gly or Asn-Ser motifs with SASA at least 30 square Angstrom
3. Oxidation sites: exposed Met (at least 30) or Trp (at least 15 square Angstrom)

Each residue is annotated with CDR membership using approximate Kabat ranges. CDR-region liabilities are prioritized as paratope-adjacent.

## Results - Trastuzumab Fab

CDR-region findings (highest priority):
- CDR-H2: clustered hydrophobic patch (Phe53, Leu54, Val58, Pro59), combined SASA over 200 square Angstrom - drives high-concentration aggregation
- CDR-H3: Pro95, Phe98 aggregation hotspots in the antigen-binding loop
- CDR-L1, CDR-L3: single hydrophobic hits (Phe27, Val93)
- CDR-L2: Asn55-Gly deamidation site (SASA = 76.8 square Angstrom), highest-priority single liability

Framework findings of note:
- Light chain Trp110: flagged for both oxidation and aggregation, consistent with published trastuzumab Trp oxidation literature

These findings are consistent with the formulation challenges that motivate Herceptin's lyophilized (rather than liquid) commercial presentation.

## Reproduce

    conda create -n biopharma -c conda-forge python=3.11 biopython freesasa pandas numpy -y
    conda activate biopharma
    wget https://files.rcsb.org/download/1N8Z.pdb
    python liability_analysis.py
    pymol highlight_liabilities.pml

## Files

- liability_analysis.py: main analysis script
- 1N8Z.pdb: input structure (trastuzumab Fab + HER2 ECD)
- liability_summary.txt: text report (CDR liabilities first)
- per_residue_liabilities.csv: full per-residue table
- highlight_liabilities.pml: auto-generated PyMOL visualization script
- trastuzumab_liabilities.png: rendered structural map of liabilities
- trastuzumab_liabilities_session.pse: saved PyMOL session

## Limitations

- Single-residue analysis; commercial tools use patch-based scores (SAP, SCM)
- Static crystal structure - does not model conformational dynamics
- No pH-dependent ionization modeling
- No prediction of viscosity, isoelectric point, or thermal stability
- Approximate Kabat numbering - full Kabat would require sequence alignment

## References

- Trastuzumab structure: Cho et al., Nature 2003 (PDB 1N8Z)
- FreeSASA: Mitternacht, F1000Research 2016
- BioPython: Cock et al., Bioinformatics 2009
- Antibody developability concepts: Jain et al., PNAS 2017
