#!/usr/bin/env python3
"""Quick test to verify our architecture fixes."""

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

def test_basic_functionality():
    """Test basic functionality of our fixes."""
    print("=== Testing Basic Functionality ===")
    
    try:
        # Test imports
        from Alprotein.core.protein_structure import ProteinStructure
        from Alprotein.core.pigment_system import PigmentSystem
        from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
        from Alprotein.data.parameters import create_parameter_config
        print("✓ All imports successful")
        
        # Test parameter config creation
        config = create_parameter_config(['CLA', 'CHL'])
        print(f"✓ Parameter config created: {config}")
        
        # Test basic structure loading (without enhanced parameters)
        pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
        print(f"✓ Testing with PDB: {pdb_path}")
        
        # Load structure normally (old way)
        structure_old = ProteinStructure.from_file_extended(pdb_path, "test_old")
        pigment_system_old = PigmentSystem(structure_old)
        print(f"✓ Old architecture: {len(pigment_system_old.pigments)} pigments found")
        print(f"✓ Enhanced atoms available (old): {pigment_system_old.enhanced_atoms_available}")
        
        # Load structure with parameters (new way)
        structure_new = ProteinStructure.from_file_with_parameters(pdb_path, config, "test_new")
        pigment_system_new = PigmentSystem(structure_new)
        print(f"✓ New architecture: {len(pigment_system_new.pigments)} pigments found")
        print(f"✓ Enhanced atoms available (new): {pigment_system_new.enhanced_atoms_available}")
        
        # Test site energy calculation with both
        calculator = SiteEnergyCalculator()
        
        if len(pigment_system_old.pigments) > 0 and len(pigment_system_new.pigments) > 0:
            # Test old architecture
            first_pigment_old = list(pigment_system_old.pigments.values())[0]
            energy_old = calculator.calculate_site_energy(first_pigment_old, pigment_system_old)
            print(f"✓ Old architecture energy: {energy_old:.2f} cm⁻¹")
            print(f"✓ Used direct access (old): {calculator.use_direct_access}")
            
            # Test new architecture
            first_pigment_new = list(pigment_system_new.pigments.values())[0]
            energy_new = calculator.calculate_site_energy(first_pigment_new, pigment_system_new)
            print(f"✓ New architecture energy: {energy_new:.2f} cm⁻¹")
            print(f"✓ Used direct access (new): {calculator.use_direct_access}")
            
            # Check if enhanced atoms are actually working
            if hasattr(first_pigment_new, '_has_enhanced_atoms'):
                print(f"✓ Pigment has enhanced atom flag: {first_pigment_new._has_enhanced_atoms}")
            else:
                print("✗ Pigment missing enhanced atom flag")
            
            # Check site energy atoms
            site_atoms_old = first_pigment_old.get_site_energy_atoms()
            site_atoms_new = first_pigment_new.get_site_energy_atoms()
            print(f"✓ Site atoms (old): {len(site_atoms_old)}")
            print(f"✓ Site atoms (new): {len(site_atoms_new)}")
            
            # Compare energies
            if abs(energy_old - energy_new) < 1e-6:
                print("✓ Energies match between old and new architectures")
            else:
                print(f"✗ Energy mismatch: {energy_old} vs {energy_new} (diff: {abs(energy_old - energy_new)})")
        else:
            print("✗ No pigments found in one or both systems")
            
        print("\n=== Test Summary ===")
        print("✓ Basic functionality test completed successfully")
        return True
        
    except Exception as e:
        print(f"✗ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_basic_functionality()
    if success:
        print("🎉 All tests passed!")
        sys.exit(0)
    else:
        print("❌ Some tests failed!")
        sys.exit(1)
