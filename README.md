# Antibody Developability Toolkit

A small, open-source Python toolkit for computational developability assessment of therapeutic antibodies. Demonstrates an end-to-end workflow: load a structure, visualize the antigen-binding interface, and identify formulation-relevant biophysical liabilities.

Built as a self-directed exploration of biopharmaceutical formulation analysis using only open-source tools (BioPython, FreeSASA, PyMOL).

## Why this exists

Therapeutic antibodies must remain stable across years of shelf life and high-concentration storage. Formulation scientists need to identify, ahead of expensive wet-lab experiments, where a candidate molecule might fail: aggregation-prone surfaces, chemically degradable residues, oxidation-sensitive sites. This toolkit demonstrates the computational side of that assessment workflow on a real, marketed therapeutic antibody.

## Case study: trastuzumab (Herceptin)

Both modules use **trastuzumab Fab bound to HER2** (PDB 1N8Z) as the demonstration molecule. Trastuzumab is one of the most successful therapeutic antibodies ever developed, treats HER2-positive breast cancer, and is sold as a lyophilized powder due to known liquid-formulation challenges - making it an excellent test case.

## Modules

### 01-structural-visualization
PyMOL-based visualization of the trastuzumab-HER2 binding interface, with CDR-H3 highlighted at the antigen-recognition site. Reproducible `.pml` script and rendered images included.

### 02-liability-analysis
Python tool combining BioPython structural parsing with FreeSASA solvent-accessibility computation to flag three classes of formulation liabilities:
- **Aggregation hotspots** - surface-exposed hydrophobic residues
- **Deamidation sites** - solvent-exposed Asn-Gly / Asn-Ser motifs
- **Oxidation sites** - exposed Met or Trp residues

Each residue is annotated with CDR membership; CDR-region findings are prioritized as paratope-adjacent.

## Key findings on trastuzumab

The analysis independently identified liabilities consistent with published trastuzumab stability literature:
- **CDR-H2 hydrophobic patch** (Phe53, Leu54, Val58, Pro59) - clustered aggregation driver
- **CDR-L2 deamidation site** (Asn55-Gly, SASA = 76.8 square Angstrom) - paratope-adjacent
- **Light-chain Trp110** - flagged for both oxidation and aggregation

These findings are consistent with the formulation challenges that motivate Herceptin's lyophilized rather than liquid commercial presentation.

## Stack

- Python 3.11
- BioPython (structural parsing)
- FreeSASA (solvent-accessible surface area)
- pandas (tabular output)
- PyMOL (visualization)
- All packages installable via conda-forge

## Limitations and what would come next

This is a demonstration toolkit, not a production developability platform. Honest gaps:
- Single-residue scoring; production tools use patch-based aggregation scores (SAP, SCM)
- Static crystal structure - no conformational dynamics
- No pH-dependent ionization, no viscosity prediction, no thermal stability modeling
- Approximate Kabat numbering

Natural extensions: integration of multiple antibodies into a feature matrix, training a developability classifier (RandomForest with SHAP) against published experimental panels (e.g. Jain et al. PNAS 2017), and adding patch-based aggregation scoring.

## Author

Chukwuma Winner Obiora - Bioinformatics scientist with a biochemistry background, exploring the application of computational pipeline engineering and machine learning to biopharmaceutical formulation development.
