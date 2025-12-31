#!/usr/bin/env python3
"""
Enhanced Example: Atom-Level CDC Analysis (User-Friendly & Publication-Ready)

This enhanced example provides an intuitive, step-by-step guide for users to:
    - Understand and run the full CDC workflow
    - Explore site energy shifts at the atom level
    - Generate publication-quality figures and tables
    - Export structured data for further analysis

Usage:
    python run_atom_level_demo.py

Output:
    - Organized logs
    - JSON and text summaries
    - Publication-ready bar plot of dominant atomic contributions
    - Table of top contributors across residues
"""

import logging
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# === Enhanced Logging === #
def setup_logging(log_dir: Path):
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "run.log"
    logger = logging.getLogger("Alprotein-AtomCDC")
    logger.setLevel(logging.DEBUG)

    if logger.hasHandlers():
        logger.handlers.clear()

    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)
    return logger


# === Visualization Helper === #
def plot_dominant_atom_bar(detailed_contributions, output_dir, cutoff=20):
    """Generate bar plot of top contributing atoms."""
    dominant_list = []
    for pigment_data in detailed_contributions.values():
        for residue_data in pigment_data.values():
            if 'dominant_atoms' in residue_data and abs(residue_data['total_contribution']) > cutoff:
                for atom, contrib in residue_data['dominant_atoms'].items():
                    dominant_list.append((atom, contrib))

    df = pd.DataFrame(dominant_list, columns=['atom', 'contribution'])
    top_atoms = df.groupby('atom')['contribution'].sum().sort_values(key=abs, ascending=False).head(10)

    plt.figure(figsize=(6, 3))
    top_atoms.plot(kind='bar', color='forestgreen')
    plt.title("Top Atom Contributions to Site Energy Shift")
    plt.ylabel("Total Contribution (cm$^{-1}$)")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_dir / "dominant_atoms_barplot.svg")
    plt.savefig(output_dir / "dominant_atoms_barplot.png", dpi=300)
    plt.close()


# === Main Workflow === #
def main():
    script_dir = Path(__file__).resolve().parent
    out_dir = script_dir / "atom_cdc_results"
    calc_dir = out_dir / "calc"
    logger = setup_logging(out_dir)

    logger.info("\n🚀 Starting Atom-Level CDC Workflow...")

    try:
        from Alprotein.core.protein_structure import ProteinStructure
        from Alprotein.core.pigment_system import PigmentSystem
        from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
        from Alprotein.utils.atom_level_cdc import (
            analyze_site_energy_contributions_with_atoms,
            find_strongest_atomic_contribution,
            analyze_backbone_vs_sidechain
        )
    except ImportError as e:
        logger.error("Import failed. Ensure Alprotein is installed correctly.\n" + str(e))
        return

    pdb_file = script_dir.parent / "Alprotein" / "pdb" / "CP24_most_occ_pH7.pdb"
    cutoff = 1000.0
    if not pdb_file.exists():
        logger.error(f"Missing PDB file: {pdb_file}")
        return

    protein = ProteinStructure.from_file_extended(str(pdb_file), 'atom_analysis')
    pigment_system = PigmentSystem(protein)
    pigment_system.add_pigments_by_residue(
        resname="CLA",
        tresp_dict_name="CLA_IPPC",
        cdc_dict_name="CLA"
    )
    pigment_system.add_pigments_by_residue(
        resname="CHL",
        tresp_dict_name="CHL_IPPC",
        cdc_dict_name="CHL"
    )
    calculator = SiteEnergyCalculator(dielectric=2.0)
    print(calculator)
    results = analyze_site_energy_contributions_with_atoms(
        pigment_system, calculator, str(calc_dir), cutoff=cutoff)
    print(results)
    logger.info("\n✅ Analysis completed. Generating outputs...")

    # Save structured output
    json_file = calc_dir / "cdc_contributions_with_atoms.json"
    txt_file = calc_dir / "cdc_summary_with_atoms.txt"
    pd.DataFrame(results).to_json(json_file, indent=2)
    with open(txt_file, 'w') as f:
        f.write("CDC Summary Table\n")
        f.write("="*60 + "\n")
        for pid, pdata in results.items():
            for rid, rdata in pdata.items():
                f.write(f"{pid}-{rid}: {rdata.get('total_contribution', 0):+.2f} cm^-1\n")

    # Generate bar plot of dominant atoms
    plot_dominant_atom_bar(results, out_dir, cutoff)
    logger.info(f"Figures saved to: {out_dir}")

    # Optional: Additional analysis
    strongest = find_strongest_atomic_contribution(results)
    if strongest:
        pid, rid, atom, value = strongest
        logger.info(f"\n🧪 Strongest Contribution: {pid}-{rid}-{atom}: {value:+.2f} cm^-1")

    bb_frac = analyze_backbone_vs_sidechain(results)
    bb_dominated = sum(1 for x in bb_frac.values() if x['backbone_fraction'] > 0.7)
    logger.info(f"\n🔍 {bb_dominated} residues show >70% backbone contribution")

    logger.info("\n🎉 Done! Explore the plots and logs for insights.")
    from pigment_site_summary import summarize_site_energy_by_pigment

    summary_df = summarize_site_energy_by_pigment(results, out_dir)
    logger.info(f"\n📊 Site energy summary saved: {out_dir / 'site_energy_shift_by_pigment.csv'}")


if __name__ == "__main__":
    main()
