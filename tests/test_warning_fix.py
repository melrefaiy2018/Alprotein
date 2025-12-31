#!/usr/bin/env python3
"""Quick test to verify the warning fix."""

import sys
import os
sys.path.insert(0, os.path.abspath('.'))

def test_warning_fix():
    """Test that warnings are eliminated."""
    print("🧪 Testing Warning Fix")
    print("=" * 30)
    
    try:
        from Alprotein.core.protein_structure import ProteinStructure
        from Alprotein.core.pigment_system import PigmentSystem
        from Alprotein.data.parameters import create_parameter_config
        
        pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
        config = create_parameter_config(['CLA', 'CHL'])
        print(f"✓ Config: {config}")
        
        print("\n--- Testing Enhanced Architecture (Should have NO warnings) ---")
        structure_new = ProteinStructure.from_file_with_parameters(pdb_path, config, "test_new")
        pigment_system_new = PigmentSystem(structure_new)
        print(f"✓ Pigments created: {len(pigment_system_new.pigments)}")
        print(f"✓ Enhanced atoms: {pigment_system_new.enhanced_atoms_available}")
        
        # Test a single pigment site energy
        if len(pigment_system_new.pigments) > 0:
            first_pigment = list(pigment_system_new.pigments.values())[0]
            site_atoms = first_pigment.get_site_energy_atoms()
            print(f"✓ Site atoms available: {len(site_atoms)}")
            print(f"✓ Has enhanced atoms flag: {getattr(first_pigment, '_has_enhanced_atoms', False)}")
        
        print("\n🎉 Test completed - check above for any warnings!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_warning_fix()
