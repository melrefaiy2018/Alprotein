"""Test script to verify migration from old to new architecture."""

from pathlib import Path
import time

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.data.parameters import create_parameter_config

def test_old_architecture(pdb_path):
    """Test using old architecture."""
    print("Testing old architecture...")
    start_time = time.time()
    
    # Load structure without parameters - this forces legacy mode
    structure = ProteinStructure.from_file_extended(pdb_path, "test")
    pigment_system = PigmentSystem(structure)
    calculator = SiteEnergyCalculator()
    
    # Calculate site energies - will use legacy path automatically
    results_old = {}
    for i, pigment in enumerate(pigment_system.pigments.values()):
        energy = calculator.calculate_site_energy(pigment, pigment_system)
        results_old[i] = energy
    
    old_time = time.time() - start_time
    print(f"Old architecture time: {old_time:.4f} seconds")
    return results_old, old_time

def test_new_architecture(pdb_path):
    """Test using new architecture."""
    print("Testing new architecture...")
    start_time = time.time()
    
    # Create parameter configuration
    temp_structure = ProteinStructure.from_file_extended(pdb_path, "temp")
    temp_system = PigmentSystem(temp_structure)
    pigment_resnames = [p.residue.resname for p in temp_system.pigments.values()]
    config = create_parameter_config(pigment_resnames)
    
    # Load structure with parameters
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
    pigment_system = PigmentSystem(structure)
    calculator = SiteEnergyCalculator()
    
    # Calculate site energies
    results_new = {}
    for i, pigment in enumerate(pigment_system.pigments.values()):
        energy = calculator.calculate_site_energy(pigment, pigment_system)
        results_new[i] = energy
    
    new_time = time.time() - start_time
    print(f"New architecture time: {new_time:.4f} seconds")
    return results_new, new_time

def compare_results(results_old, results_new, tolerance=1e-6):
    """Compare results from old and new architectures."""
    print("\nComparing results...")
    
    if len(results_old) != len(results_new):
        print(f"ERROR: Different number of results: {len(results_old)} vs {len(results_new)}")
        return False
    
    all_match = True
    for i in range(len(results_old)):
        diff = abs(results_old[i] - results_new[i])
        if diff > tolerance:
            print(f"MISMATCH at pigment {i}: {results_old[i]} vs {results_new[i]} (diff: {diff})")
            all_match = False
        else:
            print(f"MATCH at pigment {i}: {results_old[i]} vs {results_new[i]}")
    
    return all_match

if __name__ == "__main__":
    script_dir = Path(__file__).resolve().parent
    pdb_path = script_dir.parent / "Alprotein" / "pdb" / "CP24_extended_most_occ_pH8.pdb"
    
    # Test both architectures
    results_old, time_old = test_old_architecture(pdb_path)
    results_new, time_new = test_new_architecture(pdb_path)
    
    # Compare results
    match = compare_results(results_old, results_new)
    
    print(f"\nPerformance improvement: {((time_old - time_new) / time_old * 100):.1f}%")
    print(f"Results match: {match}")
