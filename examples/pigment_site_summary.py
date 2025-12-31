import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def summarize_site_energy_by_pigment(results, output_dir: Path):
    """Creates a bar plot and table for pigment-wise total site energy shift."""
    pigment_summary = []

    for pigment_id, residues in results.items():
        total_shift = sum(
            residue_data.get('total_contribution', 0.0)
            for residue_id, residue_data in residues.items()
            if isinstance(residue_data, dict)
        )
        pigment_summary.append((pigment_id, total_shift))

    df = pd.DataFrame(pigment_summary, columns=['pigment_id', 'site_energy_shift'])
    df.sort_values('site_energy_shift', key=abs, ascending=False, inplace=True)

    # Save table
    table_path = output_dir / "site_energy_shift_by_pigment.csv"
    df.to_csv(table_path, index=False)

    # Plot
    plt.figure(figsize=(6, 4))
    plt.bar(df['pigment_id'], df['site_energy_shift'], color='teal')
    plt.title("Site Energy Shift per Pigment")
    plt.ylabel("Shift (cm$^{-1}$)")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "site_energy_shift_by_pigment.svg")
    plt.savefig(output_dir / "site_energy_shift_by_pigment.png", dpi=300)
    plt.close()

    return df
