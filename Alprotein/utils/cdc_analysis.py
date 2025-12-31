"""
CDC (Chromophore-Dependent Contribution) Analysis Module

This module provides tools for analyzing and visualizing individual residue 
contributions to site energy shifts in pigment systems.
"""

import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from typing import Dict, List, Tuple, Optional, Any, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.pigment_system import PigmentSystem
    from ..calculators.site_energy_calculator import SiteEnergyCalculator


class CDCAnalyzer:
    """
    Analyzer for Chromophore-Dependent Contributions to site energy shifts.
    
    This class provides methods to calculate, analyze, and visualize how individual
    residues and atoms contribute to site energy shifts in pigment systems.
    """
    
    def __init__(self, site_energy_calculator: 'SiteEnergyCalculator'):
        """
        Initialize CDC analyzer.
        
        Args:
            site_energy_calculator: Configured site energy calculator instance
        """
        self.calculator = site_energy_calculator
    
    def calculate_detailed_contributions(self, 
                                       pigment_system: 'PigmentSystem',
                                       target_pigment_id: str) -> Dict[str, float]:
        """
        Calculate detailed residue-by-residue contributions to a pigment's site energy.
        
        Args:
            pigment_system: The pigment system
            target_pigment_id: ID of the target pigment to analyze
            
        Returns:
            Dictionary mapping residue/pigment IDs to their energy contributions
        """
        target_pigment = pigment_system.get_pigment_by_id(target_pigment_id)
        if target_pigment is None:
            raise ValueError(f"Pigment {target_pigment_id} not found in system")
        
        contributions = {}
        
        # Get vacuum energy
        vacuum_energy = target_pigment.get_vacuum_energy()
        contributions['vacuum'] = vacuum_energy
        
        # Calculate protein background contributions residue by residue
        protein_atoms = pigment_system.get_protein_background_atoms()
        residue_groups = self._group_atoms_by_residue(protein_atoms)
        
        for residue_id, atoms in residue_groups.items():
            charges = [getattr(atom, 'charge', 0.0) for atom in atoms]
            contribution = self.calculator._pigment_background_interaction(
                target_pigment, atoms, charges
            )
            if abs(contribution) > 1e-6:  # Only store non-zero contributions
                contributions[residue_id] = contribution
        
        # Calculate pigment-pigment contributions
        for other_pigment in pigment_system:
            if other_pigment.get_id() == target_pigment_id:
                continue
            
            contribution = self.calculator._pigment_pigment_interaction(
                target_pigment, other_pigment
            )
            if abs(contribution) > 1e-6:
                contributions[other_pigment.get_id()] = contribution
        
        return contributions
    
    def calculate_all_pigment_contributions(self, 
                                          pigment_system: 'PigmentSystem') -> Dict[str, Dict[str, float]]:
        """
        Calculate detailed contributions for all pigments in the system.
        
        This method now uses the improved SiteEnergyCalculator method to get
        comprehensive residue-by-residue contributions.
        
        Args:
            pigment_system: The pigment system
            
        Returns:
            Nested dictionary: {pigment_id: {contributor_id: contribution}}
        """
        print("Calculating detailed CDC contributions...")
        
        # Use the improved SiteEnergyCalculator method
        total_shifts, all_contributions = self.calculator.calculate_detailed_site_energy_contributions(pigment_system)
        
        print(f"✓ Analysis complete for {len(all_contributions)} pigments")
        return all_contributions
    
    def _group_atoms_by_residue(self, atoms: List) -> Dict[str, List]:
        """
        Group atoms by their residue identifier.
        
        Args:
            atoms: List of BioPython atoms
            
        Returns:
            Dictionary mapping residue IDs to lists of atoms
        """
        residue_groups = {}
        
        for atom in atoms:
            try:
                residue = atom.get_parent()
                chain = residue.get_parent()
                residue_id = f"{chain.id}_{residue.resname}_{residue.id[1]}"
                
                if residue_id not in residue_groups:
                    residue_groups[residue_id] = []
                residue_groups[residue_id].append(atom)
                
            except AttributeError:
                # Handle atoms without proper parent structure
                continue
        
        return residue_groups


class CDCVisualizer:
    """
    Visualization tools for CDC analysis results.
    """
    
    def __init__(self):
        """Initialize CDC visualizer."""
        pass
    
    def create_directories(self, paths: List[str]) -> None:
        """
        Create directories if they don't exist.
        
        Args:
            paths: List of directory paths to create
        """
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)
            else:
                print(f'Directory {path} already exists and may not be empty!')
    
    def save_contributions_to_file(self, file_path: str, index: str, 
                                 keys: List[str], values: List[float]) -> None:
        """
        Save contribution data to a text file.
        
        Args:
            file_path: Output file path
            index: Pigment index/name
            keys: List of contributor names
            values: List of contribution values
        """
        with open(file_path, 'w') as dict_file:
            dict_file.write(f'Residues contributing to the energy shift of {index}:\n')
            dict_file.write('=' * 60 + '\n\n')
            
            # Filter out vacuum energy - only show actual contributions
            filtered_data = [(key, value) for key, value in zip(keys, values) 
                           if key != 'vacuum']
            
            if not filtered_data:
                dict_file.write('No significant contributions above the cutoff threshold.\n')
            else:
                for key, value in filtered_data:
                    dict_file.write(f'{key}: {value:.2f} cm⁻¹\n')
    
    def plot_residue_contribution(self, contributions: Dict[str, Dict[str, float]], 
                                cutoff: float, main_path: str) -> None:
        """
        Create bar plots of residue contributions above a cutoff.
        
        Args:
            contributions: Nested dictionary of contributions
            cutoff: Minimum absolute contribution to include
            main_path: Main output directory path
        """
        residue_contribution_path = os.path.join(main_path, f'residue_contribution_{cutoff}')
        selected_contribution_path = os.path.join(main_path, f'contribution_{cutoff}')
        
        self.create_directories([residue_contribution_path, selected_contribution_path])
        
        for pigment_id, pigment_contributions in contributions.items():
            # Sort contributions by magnitude
            sorted_contributions = dict(sorted(
                pigment_contributions.items(), 
                key=lambda item: abs(item[1]), 
                reverse=True
            ))
            
            # Filter by cutoff (exclude vacuum energy from the contribution files)
            contribution_items = [(k, v) for k, v in sorted_contributions.items() 
                                if k != 'vacuum' and abs(v) > cutoff]
            
            if not contribution_items:
                continue
                
            contribution_keys = [item[0] for item in contribution_items]
            contribution_values = [item[1] for item in contribution_items]
            
            # Save to file (only actual contributions, not vacuum energy)
            file_path = os.path.join(selected_contribution_path, 
                                   f'selected_contribution_to_{pigment_id}.txt')
            self.save_contributions_to_file(file_path, pigment_id, contribution_keys, contribution_values)
            
            # Create plot (exclude vacuum energy from plot)
            plot_keys = contribution_keys  # Use the same filtered data
            plot_values = contribution_values
            
            if plot_values:
                fig, ax = plt.subplots(figsize=(12, 8), tight_layout=True)
                
                # Color code based on contribution sign
                colors = ['red' if v < 0 else 'blue' for v in plot_values]
                bars = ax.bar(range(len(plot_keys)), plot_values, color=colors, alpha=0.7)
                
                # Customize plot
                ax.set_ylabel("Energy Contribution (cm⁻¹)", fontweight='bold', fontsize=12)
                ax.set_xlabel("Contributing Residues/Pigments", fontweight='bold', fontsize=12)
                ax.set_xticks(range(len(plot_keys)))
                ax.set_xticklabels(plot_keys, rotation=45, ha='right', fontsize=10)
                ax.tick_params(axis='y', labelsize=12)
                ax.grid(axis='y', alpha=0.3)
                
                # Add value labels on bars
                for bar, value in zip(bars, plot_values):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + (5 if height > 0 else -15),
                           f'{value:.0f}', ha='center', va='bottom' if height > 0 else 'top', 
                           fontsize=9)
                
                plt.title(f"Site Energy Shift Contributions for {pigment_id}\n(cutoff = {cutoff} cm⁻¹)", 
                         fontweight='bold', fontsize=14)
                plt.tight_layout()
                plt.savefig(os.path.join(residue_contribution_path, f'pigment_{pigment_id}.png'), 
                           dpi=300, bbox_inches='tight')
                plt.close()
    
    def parse_interaction_data(self, interaction_path: str) -> List[Tuple[str, str, float]]:
        """
        Parse interaction data from a text file.
        
        Args:
            interaction_path: Path to the interaction data file
            
        Returns:
            List of (residue_pair, residue_name, interaction_value) tuples
        """
        if not os.path.exists(interaction_path):
            raise FileNotFoundError(f"The file {interaction_path} does not exist.")
        
        interaction_data = []
        with open(interaction_path, "r") as f:
            lines = f.readlines()
            
            # Skip header lines and empty lines
            for line in lines:
                line = line.strip()
                if not line or '=' in line or line.startswith('Residues contributing'):
                    continue
                if line.startswith('No significant contributions'):
                    break
                
                # Parse contribution lines
                # Expected format: "A_ASP_123: -45.67 cm⁻¹" or "Residue: A_ASP_123, Contribution: -45.67 cm⁻¹"
                if ':' in line:
                    if line.startswith('Residue:'):  # New format
                        # Format: "Residue: A_ASP_123, Contribution: -45.67 cm⁻¹"
                        parts = line.split(',')
                        if len(parts) >= 2:
                            residue_part = parts[0].replace('Residue:', '').strip()
                            contrib_part = parts[1].replace('Contribution:', '').strip()
                            try:
                                # Remove 'cm⁻¹' and convert to float
                                interaction_value = float(contrib_part.replace('cm⁻¹', '').strip())
                                interaction_data.append((residue_part, residue_part, interaction_value))
                            except (ValueError, IndexError):
                                continue
                    else:  # Simple format
                        # Format: "A_ASP_123: -45.67 cm⁻¹"
                        parts = line.split(':')
                        if len(parts) >= 2:
                            residue_name = parts[0].strip()
                            try:
                                # Remove 'cm⁻¹' and convert to float
                                value_str = parts[1].replace('cm⁻¹', '').strip()
                                interaction_value = float(value_str)
                                interaction_data.append((residue_name, residue_name, interaction_value))
                            except (ValueError, IndexError):
                                continue
        
        return interaction_data
    
    def calculate_residue_centers(self, pdb_path: str) -> Dict[str, np.ndarray]:
        """
        Calculate center of mass for each residue in a PDB file.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            Dictionary mapping residue identifiers to center positions
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        
        residue_centers = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if len(list(residue.get_atoms())) == 0:
                        continue
                    
                    center_of_mass = np.mean([atom.get_coord() for atom in residue], axis=0)
                    residue_id = f"{chain.id}_{residue.resname}_{residue.id[1]}"
                    residue_centers[residue_id] = center_of_mass
        
        return residue_centers
    
    def plot_circle_value(self, center_of_mass: np.ndarray, residue: str, 
                         interaction_value: float, color: str, circle_size: float, 
                         fontsize: int, show_values: bool, cutoff: float, 
                         shape: str = 'o') -> None:
        """
        Plot a single point with optional value label.
        
        Args:
            center_of_mass: Position coordinates
            residue: Residue identifier (full format: Chain_ResName_ResNum)
            interaction_value: Contribution value
            color: Marker color
            circle_size: Marker size
            fontsize: Font size for labels
            show_values: Whether to show value labels
            cutoff: Minimum value to show labels
            shape: Marker shape
        """
        plt.scatter(center_of_mass[0], center_of_mass[1], c=[color], s=circle_size, 
                   alpha=0.6, edgecolors="black", marker=shape, linewidth=0.5)
        
        # Only display the label if the absolute interaction value is above the cutoff
        if show_values and abs(interaction_value) > cutoff:
            # Parse residue identifier to create a readable label
            # Format: Chain_ResName_ResNum -> Chain:ResName ResNum
            try:
                parts = residue.split('_')
                if len(parts) >= 3:
                    chain_id = parts[0]
                    res_name = parts[1]
                    res_num = parts[2]
                    # Create compact label: Chain:ResName\nResNum (Value)
                    label = f"{chain_id}:{res_name}\n{res_num} ({round(interaction_value):d})"
                else:
                    # Fallback for unexpected format
                    label = f"{residue}\n{round(interaction_value):d}"
            except:
                # Fallback for any parsing errors
                label = f"{residue}\n{round(interaction_value):d}"
            
            plt.text(center_of_mass[0], center_of_mass[1], label, 
                    color="black", ha="center", va="center", fontsize=max(6, fontsize-2),
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="gray", linewidth=0.5))
    
    def plot_spatial_interaction(self, residue_centers: Dict[str, np.ndarray], 
                               interaction_data: List[Tuple[str, str, float]], 
                               interacting_pigment: str, font_size: int = 8, 
                               show_values: bool = True, main_path: str = "", 
                               cutoff: float = 0) -> None:
        """
        Create 2D spatial plot of residue contributions.
        
        Args:
            residue_centers: Dictionary of residue positions
            interaction_data: List of interaction data tuples
            interacting_pigment: Target pigment identifier
            font_size: Font size for labels
            show_values: Whether to show value labels
            main_path: Output directory path
            cutoff: Minimum contribution to visualize
        """
        residue_interactions = {entry[1]: entry[2] for entry in interaction_data}
        
        # Create larger figure to accommodate longer labels
        plt.figure(figsize=(14, 11))
        
        # Plot the target pigment in green with a distinct shape
        target_found = False
        for residue_id, center in residue_centers.items():
            if interacting_pigment in residue_id:
                self.plot_circle_value(center, interacting_pigment, 0, 'green', 
                                     800, font_size, show_values=False, 
                                     cutoff=cutoff, shape='h')
                target_found = True
                break
        
        if not target_found:
            print(f"Warning: Target pigment {interacting_pigment} not found in residue centers")
        
        # Plot contributing residues
        for residue, center in residue_centers.items():
            if interacting_pigment in residue:
                continue  # Skip the target pigment
            
            interaction_value = residue_interactions.get(residue, 0)
            if abs(interaction_value) < cutoff:
                continue
            
            # Determine color and shape based on residue type and contribution
            if 'NA' in residue:
                color, shape = 'purple', 's'
                circle_size = 300
            elif 'CL' in residue:
                color, shape = 'cyan', 's'
                circle_size = 300
            else:
                color = 'blue' if interaction_value > 0 else 'red'
                shape = 'o'
                circle_size = min(100 * abs(interaction_value), 1000)
            
            self.plot_circle_value(center, residue, interaction_value, color, 
                                 circle_size, font_size, show_values, cutoff, shape)
        
        plt.xlabel("X Coordinate (Å)", fontsize=12, fontweight='bold')
        plt.ylabel("Y Coordinate (Å)", fontsize=12, fontweight='bold')
        plt.title(f"Spatial Distribution of Energy Contributions to {interacting_pigment}\n(Chain:ResName ResNum format)", 
                 fontsize=14, fontweight='bold')
        plt.grid(alpha=0.3)
        plt.axis('equal')
        
        # Add colorbar legend with updated descriptions
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', label='Target Pigment'),
            Patch(facecolor='blue', label='Positive Contribution (Blue-shift)'),
            Patch(facecolor='red', label='Negative Contribution (Red-shift)'),
            Patch(facecolor='purple', label='Na⁺ Ions'),
            Patch(facecolor='cyan', label='Cl⁻ Ions')
        ]
        plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1), 
                  fontsize=10, framealpha=0.9)
        
        plt.tight_layout()
        
        if main_path:
            save_path_fig = os.path.join(main_path, '2d_interaction_fig')
            save_path_data = os.path.join(main_path, '2d_interaction_data')
            
            self.create_directories([save_path_fig, save_path_data])
            plt.savefig(os.path.join(save_path_fig, f'spatial_interaction_{interacting_pigment}.png'), 
                       dpi=600, bbox_inches='tight')
            
            # Save interaction data
            with open(os.path.join(save_path_data, f'residue_interactions_{interacting_pigment}.txt'), 'w') as txt_file:
                txt_file.write(f'Residues contributing to the energy shift of {interacting_pigment}:\n')
                txt_file.write('=' * 60 + '\n\n')
                
                # Filter out vacuum energy and sort by magnitude
                filtered_interactions = {k: v for k, v in residue_interactions.items() 
                                       if k != 'vacuum'}
                sorted_interactions = sorted(filtered_interactions.items(), 
                                           key=lambda x: abs(x[1]), reverse=True)
                
                if not sorted_interactions:
                    txt_file.write('No significant contributions above the cutoff threshold.\n')
                else:
                    for residue, interaction_value in sorted_interactions:
                        if abs(interaction_value) > cutoff:
                            txt_file.write(f"Residue: {residue}, Contribution: {interaction_value:.2f} cm⁻¹\n")
        
        plt.close()


class CDCExporter:
    """
    Export and summary tools for CDC analysis.
    """
    
    def __init__(self):
        """Initialize CDC exporter."""
        pass
    
    def create_summary_report(self, contributions: Dict[str, Dict[str, float]], 
                             cutoff: float, output_path: str) -> None:
        """
        Create a comprehensive summary report of all CDC contributions.
        
        Args:
            contributions: Nested dictionary of contributions
            cutoff: Minimum contribution threshold
            output_path: Output file path
        """
        with open(output_path, 'w') as f:
            f.write("CDC Analysis Summary Report\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Analysis performed with cutoff: {cutoff} cm⁻¹\n")
            f.write(f"Total pigments analyzed: {len(contributions)}\n\n")
            
            for pigment_id, pigment_contribs in contributions.items():
                f.write(f"\n{pigment_id}:\n")
                f.write("-" * (len(pigment_id) + 1) + "\n")
                
                # Separate vacuum energy and contributions
                vacuum = pigment_contribs.get('vacuum', 0.0)
                actual_contribs = {k: v for k, v in pigment_contribs.items() 
                                 if k != 'vacuum' and abs(v) > cutoff}
                
                f.write(f"  Vacuum Energy: {vacuum:.1f} cm⁻¹\n")
                
                if actual_contribs:
                    total_shift = sum(actual_contribs.values())
                    f.write(f"  Total Energy Shift: {total_shift:+.1f} cm⁻¹\n")
                    f.write(f"  Final Site Energy: {vacuum + total_shift:.1f} cm⁻¹\n\n")
                    
                    # Sort contributions by magnitude
                    sorted_contribs = sorted(actual_contribs.items(), 
                                           key=lambda x: abs(x[1]), reverse=True)
                    
                    f.write("  Major Contributors:\n")
                    for contributor, value in sorted_contribs[:10]:  # Top 10
                        f.write(f"    {contributor}: {value:+.1f} cm⁻¹\n")
                    
                    if len(sorted_contribs) > 10:
                        f.write(f"    ... and {len(sorted_contribs) - 10} more\n")
                else:
                    f.write("  No significant contributions above cutoff\n")
                
                f.write("\n")
    
    def create_summary_report_with_atoms(self, contributions_with_atoms: Dict[str, Dict[str, Any]], 
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
                                  if k != 'vacuum' and 'total_contribution' in v 
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
    
    def export_to_json(self, contributions: Dict[str, Dict[str, float]], 
                      output_path: str) -> None:
        """
        Export contributions to JSON format.
        
        Args:
            contributions: Contribution data
            output_path: Output JSON file path
        """
        import json
        
        # Convert numpy arrays to lists if present
        export_data = {}
        for pigment_id, contribs in contributions.items():
            export_data[pigment_id] = {k: float(v) for k, v in contribs.items()}
        
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=2)


def analyze_site_energy_contributions(pigment_system: 'PigmentSystem',
                                    site_energy_calculator: 'SiteEnergyCalculator',
                                    pdb_path: str,
                                    output_dir: str,
                                    cutoff: float = 20.0,
                                    font_size: int = 8,
                                    show_values: bool = True) -> Dict[str, Dict[str, float]]:
    """
    Complete CDC analysis workflow.
    
    Args:
        pigment_system: The pigment system to analyze
        site_energy_calculator: Configured site energy calculator
        pdb_path: Path to the PDB file for spatial visualization
        output_dir: Directory to save results
        cutoff: Minimum contribution threshold (cm⁻¹)
        font_size: Font size for plots
        show_values: Whether to show values on spatial plots
        
    Returns:
        Dictionary of all calculated contributions
    """
    print("Starting comprehensive CDC analysis...")
    
    # Initialize analyzers
    analyzer = CDCAnalyzer(site_energy_calculator)
    visualizer = CDCVisualizer()
    exporter = CDCExporter()
    
    # Create main output directory
    main_path = os.path.join(output_dir, 'CDC_analysis')
    visualizer.create_directories([main_path])
    
    # Calculate all contributions
    all_contributions = analyzer.calculate_all_pigment_contributions(pigment_system)
    
    # Create visualizations
    print("Creating contribution plots...")
    visualizer.plot_residue_contribution(all_contributions, cutoff, main_path)
    
    # Create spatial plots
    if os.path.exists(pdb_path):
        print("Creating spatial interaction plots...")
        residue_centers = visualizer.calculate_residue_centers(pdb_path)
        
        for pigment_id in all_contributions.keys():
            contribution_file = os.path.join(main_path, f'contribution_{cutoff}', 
                                           f'selected_contribution_to_{pigment_id}.txt')
            
            if os.path.exists(contribution_file):
                try:
                    interaction_data = visualizer.parse_interaction_data(contribution_file)
                    visualizer.plot_spatial_interaction(
                        residue_centers, interaction_data, pigment_id, 
                        font_size, show_values, main_path, cutoff
                    )
                except Exception as e:
                    print(f"Warning: Could not create spatial plot for {pigment_id}: {e}")
    else:
        print(f"Warning: PDB file {pdb_path} not found. Skipping spatial plots.")
    
    # Create summary reports
    print("Generating summary reports...")
    summary_path = os.path.join(main_path, 'cdc_summary_report.txt')
    exporter.create_summary_report(all_contributions, cutoff, summary_path)
    
    json_path = os.path.join(main_path, 'cdc_contributions.json')
    exporter.export_to_json(all_contributions, json_path)
    
    print(f"CDC analysis complete. Results saved to: {main_path}")
    
    return all_contributions
