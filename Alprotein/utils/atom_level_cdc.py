#!/usr/bin/env python3
"""
Analyze site energy contributions with atom-level breakdown.

This module provides enhanced CDC analysis that shows which specific atoms
within each residue contribute most to the site energy shifts.
"""

import os
import json
from typing import Dict, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ..calculators.site_energy_calculator import SiteEnergyCalculator
    from ..core.pigment_system import PigmentSystem


def analyze_site_energy_contributions_with_atoms(
    pigment_system: 'PigmentSystem',
    site_energy_calculator: 'SiteEnergyCalculator',
    output_dir: str,
    cutoff: float = 20.0
) -> Dict[str, Dict[str, Any]]:
    """
    Perform comprehensive CDC analysis with atom-level breakdown.
    
    Args:
        pigment_system: The pigment system to analyze
        site_energy_calculator: Configured site energy calculator
        output_dir: Directory to save results
        cutoff: Minimum contribution threshold (cm⁻¹)
        
    Returns:
        Dictionary with detailed atom-level contributions
    """
    print("Starting comprehensive CDC analysis with atom-level breakdown...")
    
    # Calculate detailed contributions with atom breakdown
    total_shifts, detailed_contributions = site_energy_calculator.calculate_detailed_site_energy_contributions_with_atoms(pigment_system)
    
    # Create output directory
    main_path = os.path.join(output_dir, 'CDC_analysis_with_atoms')
    if not os.path.exists(main_path):
        os.makedirs(main_path)
    
    # Generate comprehensive report with atom details
    summary_path = os.path.join(main_path, 'cdc_summary_with_atoms.txt')
    create_summary_report_with_atoms(detailed_contributions, cutoff, summary_path)
    
    # Export to JSON for further analysis
    json_path = os.path.join(main_path, 'cdc_contributions_with_atoms.json')
    export_atom_data_to_json(detailed_contributions, json_path)
    
    # Display summary
    print(f"\n📊 Atom-Level CDC Analysis Results:")
    print("-" * 50)
    
    for pigment_id, pigment_data in detailed_contributions.items():
        vacuum_energy = pigment_data.get('vacuum', {}).get('total_contribution', 0.0)
        
        # Count residues with significant contributions
        significant_residues = {k: v for k, v in pigment_data.items() 
                              if k != 'vacuum' and isinstance(v, dict) and 'total_contribution' in v 
                              and abs(v['total_contribution']) > cutoff}
        
        if significant_residues:
            total_shift = sum(v['total_contribution'] for v in significant_residues.values())
            
            print(f"\n{pigment_id}:")
            print(f"  Vacuum Energy: {vacuum_energy:.1f} cm⁻¹")
            print(f"  Total Shift: {total_shift:+.1f} cm⁻¹")
            print(f"  Final Site Energy: {vacuum_energy + total_shift:.1f} cm⁻¹")
            print(f"  Contributing residues: {len(significant_residues)}")
            
            # Show top residues with their dominant atoms
            sorted_residues = sorted(significant_residues.items(), 
                                   key=lambda x: abs(x[1]['total_contribution']), reverse=True)
            
            print(f"  Top contributors:")
            for i, (residue_id, residue_data) in enumerate(sorted_residues[:3]):
                total_contrib = residue_data['total_contribution']
                dominant_atoms = residue_data.get('dominant_atoms', {})
                
                print(f"    {residue_id}: {total_contrib:+.2f} cm⁻¹")
                
                if dominant_atoms:
                    top_atom = list(dominant_atoms.items())[0]  # Most dominant atom
                    atom_name, atom_contrib = top_atom
                    print(f"      → Dominant atom: {atom_name} ({atom_contrib:+.2f} cm⁻¹)")
                    
                    if len(dominant_atoms) > 1:
                        print(f"      → {len(dominant_atoms)} total significant atoms")
    
    print(f"\n✅ Atom-level analysis complete!")
    print(f"📁 Results saved to: {main_path}")
    print(f"📄 Detailed report: {summary_path}")
    
    return detailed_contributions


def create_summary_report_with_atoms(contributions_with_atoms: Dict[str, Dict[str, Any]], 
                                    cutoff: float, output_path: str) -> None:
    """
    Create a comprehensive summary report with atom-level breakdown.
    
    Args:
        contributions_with_atoms: Nested dictionary with atom-level contributions
        cutoff: Minimum contribution threshold
        output_path: Output file path
    """
    with open(output_path, 'w') as f:
        f.write("CDC Analysis Summary Report with Atom-Level Breakdown\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Analysis performed with cutoff: {cutoff} cm⁻¹\n")
        f.write(f"Total pigments analyzed: {len(contributions_with_atoms)}\n\n")
        
        for pigment_id, pigment_data in contributions_with_atoms.items():
            f.write(f"\n{pigment_id}:\n")
            f.write("-" * (len(pigment_id) + 1) + "\n")
            
            # Get vacuum energy
            vacuum_data = pigment_data.get('vacuum', {})
            vacuum_energy = vacuum_data.get('total_contribution', 0.0)
            f.write(f"  Vacuum Energy: {vacuum_energy:.1f} cm⁻¹\n")
            
            # Process residue contributions
            residue_contribs = {k: v for k, v in pigment_data.items() 
                              if k != 'vacuum' and isinstance(v, dict) and 'total_contribution' in v 
                              and abs(v['total_contribution']) > cutoff}
            
            if residue_contribs:
                total_shift = sum(v['total_contribution'] for v in residue_contribs.values())
                f.write(f"  Total Energy Shift: {total_shift:+.1f} cm⁻¹\n")
                f.write(f"  Final Site Energy: {vacuum_energy + total_shift:.1f} cm⁻¹\n\n")
                
                # Sort residues by total contribution magnitude
                sorted_residues = sorted(residue_contribs.items(), 
                                       key=lambda x: abs(x[1]['total_contribution']), reverse=True)
                
                f.write("  Major Contributing Residues:\n")
                for residue_id, residue_data in sorted_residues[:10]:  # Top 10
                    total_contrib = residue_data['total_contribution']
                    f.write(f"\n    {residue_id}: {total_contrib:+.2f} cm⁻¹\n")
                    
                    # Show dominant atoms within this residue
                    dominant_atoms = residue_data.get('dominant_atoms', {})
                    if dominant_atoms:
                        f.write(f"      Dominant atoms:\n")
                        for atom_name, atom_contrib in list(dominant_atoms.items())[:5]:  # Top 5 atoms
                            f.write(f"        {atom_name}: {atom_contrib:+.2f} cm⁻¹\n")
                        
                        if len(dominant_atoms) > 5:
                            f.write(f"        ... and {len(dominant_atoms) - 5} more atoms\n")
                    else:
                        f.write(f"      (No dominant atoms above threshold)\n")
                
                if len(sorted_residues) > 10:
                    f.write(f"\n    ... and {len(sorted_residues) - 10} more residues\n")
            else:
                f.write("  No significant contributions above cutoff\n")
            
            f.write("\n")


def export_atom_data_to_json(contributions_with_atoms: Dict[str, Dict[str, Any]], output_path: str) -> None:
    """Export atom-level contributions to JSON format."""
    
    # Convert to JSON-serializable format
    export_data = {}
    for pigment_id, pigment_data in contributions_with_atoms.items():
        export_data[pigment_id] = {}
        
        for residue_id, residue_data in pigment_data.items():
            if isinstance(residue_data, dict):
                # Convert all float values to ensure JSON compatibility
                serializable_data = {}
                for key, value in residue_data.items():
                    if isinstance(value, dict):
                        serializable_data[key] = {k: float(v) if isinstance(v, (int, float)) else v 
                                                for k, v in value.items()}
                    else:
                        serializable_data[key] = float(value) if isinstance(value, (int, float)) else value
                
                export_data[pigment_id][residue_id] = serializable_data
            else:
                export_data[pigment_id][residue_id] = float(residue_data) if isinstance(residue_data, (int, float)) else residue_data
    
    with open(output_path, 'w') as f:
        json.dump(export_data, f, indent=2)


def display_atom_breakdown(pigment_id: str, residue_id: str, contributions_with_atoms: Dict[str, Dict[str, Any]]) -> None:
    """Display detailed atom breakdown for a specific residue."""
    
    pigment_data = contributions_with_atoms.get(pigment_id, {})
    residue_data = pigment_data.get(residue_id, {})
    
    if not residue_data or 'atom_breakdown' not in residue_data:
        print(f"No atom breakdown available for {pigment_id} - {residue_id}")
        return
    
    total_contrib = residue_data['total_contribution']
    atom_breakdown = residue_data['atom_breakdown']
    dominant_atoms = residue_data.get('dominant_atoms', {})
    
    print(f"\n🔬 Atom-Level Breakdown for {pigment_id} - {residue_id}")
    print(f"Total Residue Contribution: {total_contrib:+.2f} cm⁻¹")
    print("-" * 50)
    
    # Sort atoms by contribution magnitude
    sorted_atoms = sorted(atom_breakdown.items(), key=lambda x: abs(x[1]), reverse=True)
    
    print("All atoms:")
    for atom_name, atom_contrib in sorted_atoms:
        is_dominant = atom_name in dominant_atoms
        marker = "★" if is_dominant else " "
        print(f"  {marker} {atom_name}: {atom_contrib:+.2f} cm⁻¹")
    
    if dominant_atoms:
        print(f"\n★ = Dominant atoms (> 10% of total or > 5 cm⁻¹)")


def find_residues_with_dominant_atom_type(contributions: Dict[str, Dict[str, Any]], atom_type: str = 'O'):
    """Find residues where a specific atom type dominates."""
    results = []
    for pigment_id, pigment_data in contributions.items():
        for residue_id, residue_data in pigment_data.items():
            if residue_id == 'vacuum' or not isinstance(residue_data, dict):
                continue
            
            dominant_atoms = residue_data.get('dominant_atoms', {})
            for atom_name, atom_contrib in dominant_atoms.items():
                if atom_type in atom_name:  # e.g., 'O1', 'OE1', etc.
                    results.append((pigment_id, residue_id, atom_name, atom_contrib))
    
    return sorted(results, key=lambda x: abs(x[3]), reverse=True)


def find_strongest_atomic_contribution(contributions: Dict[str, Dict[str, Any]]):
    """Get the most contributing atom across all residues."""
    strongest = None
    max_contrib = 0
    
    for pigment_id, pigment_data in contributions.items():
        for residue_id, residue_data in pigment_data.items():
            if residue_id == 'vacuum' or not isinstance(residue_data, dict):
                continue
                
            atom_breakdown = residue_data.get('atom_breakdown', {})
            for atom_name, atom_contrib in atom_breakdown.items():
                if abs(atom_contrib) > max_contrib:
                    max_contrib = abs(atom_contrib)
                    strongest = (pigment_id, residue_id, atom_name, atom_contrib)
    
    return strongest


def analyze_backbone_vs_sidechain(contributions: Dict[str, Dict[str, Any]]):
    """Analyze backbone vs side chain contributions."""
    backbone_atoms = ['N', 'CA', 'C', 'O']
    results = {}
    
    for pigment_id, pigment_data in contributions.items():
        for residue_id, residue_data in pigment_data.items():
            if residue_id == 'vacuum' or not isinstance(residue_data, dict):
                continue
                
            atom_breakdown = residue_data.get('atom_breakdown', {})
            
            backbone_total = sum(v for k, v in atom_breakdown.items() if k in backbone_atoms)
            sidechain_total = sum(v for k, v in atom_breakdown.items() if k not in backbone_atoms)
            
            total = backbone_total + sidechain_total
            results[f"{pigment_id}_{residue_id}"] = {
                'backbone': backbone_total,
                'sidechain': sidechain_total,
                'backbone_fraction': backbone_total / total if total != 0 else 0
            }
    
    return results


if __name__ == "__main__":
    print("This module provides atom-level CDC analysis functions.")
    print("Use analyze_site_energy_contributions_with_atoms() for comprehensive analysis.")
    print("\nAvailable functions:")
    print("  • analyze_site_energy_contributions_with_atoms() - Main analysis function")
    print("  • display_atom_breakdown() - Show detailed atom contributions")
    print("  • find_residues_with_dominant_atom_type() - Find atom type patterns")
    print("  • find_strongest_atomic_contribution() - Find most contributing atom")
    print("  • analyze_backbone_vs_sidechain() - Compare backbone/sidechain contributions")
