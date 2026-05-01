"""
Biopharmaceutical Liability Analysis — Trastuzumab Fab (PDB 1N8Z)

Computes three classes of formulation liabilities from a protein structure:
  1. Aggregation hotspots — surface-exposed hydrophobic residues
  2. Deamidation sites    — N-G or N-S motifs (solvent-exposed)
  3. Oxidation sites      — exposed Met or Trp (lower threshold for Trp)

Annotates each residue with CDR membership using approximate Kabat ranges.

Outputs:
  - per_residue_liabilities.csv     : full table including CDR annotation
  - liability_summary.txt           : flagged residues report (CDRs first)
  - highlight_liabilities.pml       : PyMOL script to visualize hits on structure
"""

import freesasa
import pandas as pd
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
import warnings
warnings.filterwarnings('ignore')

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

PDB_FILE = "1N8Z.pdb"
ANTIBODY_CHAINS = ["A", "B"]  # A = heavy, B = light, C = HER2 (ignored)

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}

# CDR ranges — approximate Kabat numbering for an IgG1 Fab
# These are used as rough flags; real Kabat would need full numbering.
CDR_RANGES = {
    "A": {  # heavy chain
        "CDR-H1": (26, 35),
        "CDR-H2": (50, 65),
        "CDR-H3": (95, 102),
    },
    "B": {  # light chain
        "CDR-L1": (24, 34),
        "CDR-L2": (50, 56),
        "CDR-L3": (89, 97),
    },
}

# SASA thresholds (Å²)
SASA_HYDROPHOBIC_THRESHOLD = 20.0
SASA_DEAMIDATION_THRESHOLD = 30.0
SASA_MET_THRESHOLD = 30.0
SASA_TRP_THRESHOLD = 15.0   # lower — Trp is large; partial exposure still oxidation-relevant

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def annotate_cdr(chain_id, resnum):
    """Return CDR name if residue falls in a CDR range, else 'framework'."""
    if chain_id not in CDR_RANGES:
        return "non-antibody"
    for cdr_name, (start, end) in CDR_RANGES[chain_id].items():
        if start <= resnum <= end:
            return cdr_name
    return "framework"

# -----------------------------------------------------------------------------
# Step 1: SASA via freesasa
# -----------------------------------------------------------------------------

print(f"Loading structure: {PDB_FILE}")
structure_freesasa = freesasa.Structure(PDB_FILE)
result = freesasa.calc(structure_freesasa)
residue_sasa = result.residueAreas()

# -----------------------------------------------------------------------------
# Step 2: Parse with BioPython, build records
# -----------------------------------------------------------------------------

print("Parsing structure with BioPython")
parser = PDBParser(QUIET=True)
structure_bio = parser.get_structure("trastuzumab", PDB_FILE)

records = []
for model in structure_bio:
    for chain in model:
        chain_id = chain.id
        if chain_id not in ANTIBODY_CHAINS:
            continue

        residues = [r for r in chain if is_aa(r, standard=True)]

        for i, residue in enumerate(residues):
            res_name = residue.get_resname()
            res_num = residue.get_id()[1]
            one_letter = seq1(res_name)

            sasa_total = None
            if chain_id in residue_sasa:
                key = str(res_num)
                if key in residue_sasa[chain_id]:
                    sasa_total = residue_sasa[chain_id][key].total

            next_res_name = residues[i + 1].get_resname() if i + 1 < len(residues) else None

            records.append({
                "chain": chain_id,
                "resnum": res_num,
                "resname": res_name,
                "aa": one_letter,
                "sasa_total": sasa_total,
                "next_resname": next_res_name,
                "region": annotate_cdr(chain_id, res_num),
            })

df = pd.DataFrame(records)
print(f"Parsed {len(df)} antibody residues across chains {ANTIBODY_CHAINS}")

# -----------------------------------------------------------------------------
# Step 3: Flag liability classes
# -----------------------------------------------------------------------------

df["is_aggregation_hotspot"] = (
    df["resname"].isin(HYDROPHOBIC)
    & (df["sasa_total"] >= SASA_HYDROPHOBIC_THRESHOLD)
)

df["is_deamidation_site"] = (
    (df["resname"] == "ASN")
    & (df["next_resname"].isin({"GLY", "SER"}))
    & (df["sasa_total"] >= SASA_DEAMIDATION_THRESHOLD)
)

# Oxidation with residue-specific thresholds
def is_oxidation(row):
    if row["resname"] == "MET" and (row["sasa_total"] or 0) >= SASA_MET_THRESHOLD:
        return True
    if row["resname"] == "TRP" and (row["sasa_total"] or 0) >= SASA_TRP_THRESHOLD:
        return True
    return False

df["is_oxidation_site"] = df.apply(is_oxidation, axis=1)

df["any_liability"] = (
    df["is_aggregation_hotspot"]
    | df["is_deamidation_site"]
    | df["is_oxidation_site"]
)

# -----------------------------------------------------------------------------
# Step 4: Save CSV
# -----------------------------------------------------------------------------

df.to_csv("per_residue_liabilities.csv", index=False)
print("Wrote per_residue_liabilities.csv")

# -----------------------------------------------------------------------------
# Step 5: Summary report — CDR liabilities prioritized
# -----------------------------------------------------------------------------

cdr_liabilities = df[(df["region"].str.startswith("CDR")) & df["any_liability"]]

agg = df[df["is_aggregation_hotspot"]]
deam = df[df["is_deamidation_site"]]
oxid = df[df["is_oxidation_site"]]

with open("liability_summary.txt", "w") as f:
    f.write("Trastuzumab Fab — Biophysical Liability Analysis\n")
    f.write(f"PDB: 1N8Z  |  Chains: {', '.join(ANTIBODY_CHAINS)}\n")
    f.write("=" * 60 + "\n\n")

    f.write("CDR LIABILITIES (highest concern — paratope-adjacent)\n")
    f.write("-" * 60 + "\n")
    if len(cdr_liabilities) == 0:
        f.write("  None.\n")
    else:
        for _, r in cdr_liabilities.iterrows():
            flags = []
            if r["is_aggregation_hotspot"]: flags.append("AGG")
            if r["is_deamidation_site"]:    flags.append("DEAM")
            if r["is_oxidation_site"]:      flags.append("OX")
            flag_str = ",".join(flags)
            f.write(f"  {r['region']:<8} Chain {r['chain']}  {r['resname']}{r['resnum']:<4}  "
                    f"SASA={r['sasa_total']:>6.1f}  [{flag_str}]\n")
    f.write("\n")

    f.write(f"AGGREGATION HOTSPOTS — surface-exposed hydrophobic ({len(agg)})\n")
    f.write("-" * 60 + "\n")
    for _, r in agg.iterrows():
        f.write(f"  Chain {r['chain']}  {r['resname']}{r['resnum']:<4}  "
                f"SASA={r['sasa_total']:>6.1f}  ({r['region']})\n")
    f.write("\n")

    f.write(f"DEAMIDATION SITES — exposed N-G or N-S ({len(deam)})\n")
    f.write("-" * 60 + "\n")
    for _, r in deam.iterrows():
        f.write(f"  Chain {r['chain']}  ASN{r['resnum']:<4} -> {r['next_resname']}  "
                f"SASA={r['sasa_total']:>6.1f}  ({r['region']})\n")
    f.write("\n")

    f.write(f"OXIDATION SITES — exposed Met or Trp ({len(oxid)})\n")
    f.write("-" * 60 + "\n")
    for _, r in oxid.iterrows():
        f.write(f"  Chain {r['chain']}  {r['resname']}{r['resnum']:<4}  "
                f"SASA={r['sasa_total']:>6.1f}  ({r['region']})\n")

print("Wrote liability_summary.txt")

# -----------------------------------------------------------------------------
# Step 6: Generate PyMOL visualization script
# -----------------------------------------------------------------------------

with open("highlight_liabilities.pml", "w") as f:
    f.write("# Auto-generated PyMOL script\n")
    f.write("# Visualizes biophysical liabilities on trastuzumab Fab\n\n")

    f.write("fetch 1n8z\n")
    f.write("remove resn HOH\n")
    f.write("bg_color white\n")
    f.write("color grey80, chain A or chain B\n")
    f.write("color palegreen, chain C\n\n")

    # Aggregation = orange spheres on CA
    agg_sel = " or ".join(
        f"(chain {r['chain']} and resi {r['resnum']})"
        for _, r in agg.iterrows()
    )
    if agg_sel:
        f.write(f"select aggregation, {agg_sel}\n")
        f.write("color orange, aggregation\n")
        f.write("show sticks, aggregation and not name N+C+O\n\n")

    # Deamidation = red sticks
    deam_sel = " or ".join(
        f"(chain {r['chain']} and resi {r['resnum']})"
        for _, r in deam.iterrows()
    )
    if deam_sel:
        f.write(f"select deamidation, {deam_sel}\n")
        f.write("color red, deamidation\n")
        f.write("show sticks, deamidation\n\n")

    # Oxidation = magenta sticks
    oxid_sel = " or ".join(
        f"(chain {r['chain']} and resi {r['resnum']})"
        for _, r in oxid.iterrows()
    )
    if oxid_sel:
        f.write(f"select oxidation, {oxid_sel}\n")
        f.write("color magenta, oxidation\n")
        f.write("show sticks, oxidation\n\n")

    f.write("zoom chain A or chain B\n")
    f.write("ray 1600, 1200\n")
    f.write("png trastuzumab_liabilities.png\n")
    f.write("save trastuzumab_liabilities_session.pse\n")

print("Wrote highlight_liabilities.pml")
print("\nDone.")
