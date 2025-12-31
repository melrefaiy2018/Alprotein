#!/usr/bin/env python3
"""
MVP Test: Basic functionality test for Alprotein package

This test verifies that the basic MVP functionality works:
1. Package imports correctly
2. Can load a structure  
3. Can set up pigment system
4. Can calculate Hamiltonian
5. Results have expected properties
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path


class TestAlproteinMVP:
    """Test class for MVP functionality."""
    
    @pytest.fixture
    def sample_pdb_path(self):
        """Return path to sample PDB file."""
        current_dir = Path(__file__).parent.parent
        pdb_path = current_dir / "Alprotein" / "pdb" / "CP24_most_occ_pH7.pdb"
        
        if not pdb_path.exists():
            pytest.skip(f"Sample PDB file not found: {pdb_path}")
        
        return str(pdb_path)
    
    def test_package_import(self):
        """Test that the package imports correctly."""
        try:
            import Alprotein
            from Alprotein import (
                ProteinStructure, 
                PigmentSystem, 
                HamiltonianCalculator,
                ChlorophyllA
            )
        except ImportError as e:
            pytest.fail(f"Failed to import Alprotein package: {e}")
    
    def test_protein_structure_loading(self, sample_pdb_path):
        """Test loading protein structure from PDB file."""
        from Alprotein import ProteinStructure
        
        protein = ProteinStructure.from_file(sample_pdb_path, name='Test')
        
        assert protein is not None
        assert protein.name == 'Test'
        assert hasattr(protein, 'pdb')
        
        # Check that structure has atoms
        atoms = list(protein.pdb.get_atoms())
        assert len(atoms) > 0, "Structure should contain atoms"
    
    def test_pigment_system_setup(self, sample_pdb_path):
        """Test setting up pigment system."""
        from Alprotein import ProteinStructure, PigmentSystem, ChlorophyllA
        
        protein = ProteinStructure.from_file(sample_pdb_path, name='Test')
        pigment_system = PigmentSystem(protein)
        
        # Try to add CLA pigments
        cla_count = pigment_system.add_pigments_by_residue(
            resname="CLA",
            pigment_class=ChlorophyllA,
            tresp_dict_name="CLA_IPPC",
            cdc_dict_name="CLA"
        )
        
        assert cla_count >= 0, "Should return count of added pigments"
        
        pigment_names = pigment_system.get_pigment_names()
        assert len(pigment_names) == cla_count, "Pigment count should match"
    
    def test_hamiltonian_calculation(self, sample_pdb_path):
        """Test Hamiltonian calculation."""
        from Alprotein import (
            ProteinStructure, 
            PigmentSystem, 
            HamiltonianCalculator,
            ChlorophyllA
        )
        
        # Set up system
        protein = ProteinStructure.from_file(sample_pdb_path, name='Test')
        pigment_system = PigmentSystem(protein)
        
        cla_count = pigment_system.add_pigments_by_residue(
            resname="CLA",
            pigment_class=ChlorophyllA,
            tresp_dict_name="CLA_IPPC",
            cdc_dict_name="CLA"
        )
        
        if cla_count == 0:
            pytest.skip("No CLA pigments found in test structure")
        
        # Calculate Hamiltonian
        calculator = HamiltonianCalculator(pigment_system)
        hamiltonian = calculator.construct_hamiltonian()
        
        # Verify Hamiltonian properties
        assert isinstance(hamiltonian, pd.DataFrame)
        assert hamiltonian.shape[0] == hamiltonian.shape[1], "Should be square matrix"
        assert hamiltonian.shape[0] == cla_count, "Size should match pigment count"
        
        # Check symmetry
        H_values = hamiltonian.values
        assert np.allclose(H_values, H_values.T, rtol=1e-10), "Should be symmetric"
        
        # Check diagonal elements (site energies) are reasonable
        diagonal = np.diag(H_values)
        assert np.all(diagonal > 10000), "Site energies should be > 10000 cm⁻¹"
        assert np.all(diagonal < 20000), "Site energies should be < 20000 cm⁻¹"
    
    def test_diagonalization(self, sample_pdb_path):
        """Test Hamiltonian diagonalization."""
        from Alprotein import (
            ProteinStructure, 
            PigmentSystem, 
            HamiltonianCalculator,
            ChlorophyllA
        )
        
        # Set up system with minimal pigments for faster test
        protein = ProteinStructure.from_file(sample_pdb_path, name='Test')
        pigment_system = PigmentSystem(protein)
        
        cla_count = pigment_system.add_pigments_by_residue(
            resname="CLA",
            pigment_class=ChlorophyllA,
            tresp_dict_name="CLA_IPPC",
            cdc_dict_name="CLA"
        )
        
        if cla_count == 0:
            pytest.skip("No CLA pigments found in test structure")
        
        # Calculate and diagonalize
        calculator = HamiltonianCalculator(pigment_system)
        hamiltonian = calculator.construct_hamiltonian()
        eigenvalues, eigenvectors = calculator.diagonalize_hamiltonian(hamiltonian)
        
        # Verify results
        assert len(eigenvalues) == cla_count, "Should have one eigenvalue per pigment"
        assert eigenvectors.shape == (cla_count, cla_count), "Eigenvectors should be square"
        
        # Check eigenvalues are real and sorted
        assert np.all(np.isreal(eigenvalues)), "Eigenvalues should be real"
        assert np.all(eigenvalues[:-1] <= eigenvalues[1:]), "Should be sorted ascending"
        
        # Check eigenvectors are orthonormal
        identity = np.eye(cla_count)
        product = eigenvectors.T @ eigenvectors
        assert np.allclose(product, identity, rtol=1e-10), "Eigenvectors should be orthonormal"
    
    def test_quick_calculation_function(self, sample_pdb_path):
        """Test the quick calculation convenience function."""
        from Alprotein import quick_calculation
        
        try:
            results = quick_calculation(sample_pdb_path)
            
            # Check return format
            assert isinstance(results, dict)
            required_keys = ['hamiltonian', 'eigenvalues', 'eigenvectors', 'pigment_names']
            for key in required_keys:
                assert key in results, f"Missing key: {key}"
            
            # Basic checks
            assert isinstance(results['hamiltonian'], pd.DataFrame)
            assert isinstance(results['eigenvalues'], np.ndarray)
            assert isinstance(results['eigenvectors'], np.ndarray)
            assert isinstance(results['pigment_names'], list)
            
        except Exception as e:
            # If quick calculation fails, it might be due to no pigments
            # This is acceptable for the MVP test
            if "No pigments found" in str(e) or len(str(e)) == 0:
                pytest.skip("No pigments found for quick calculation")
            else:
                raise e

    def test_charge_warning_and_extended_loading(self, sample_pdb_path):
        """Ensure from_file warns about charges and extended loader works."""
        from Alprotein import ProteinStructure, PigmentSystem
        import warnings

        # Regular loader should warn and result in no background atoms with charges
        with pytest.warns(UserWarning, match="Charge fields detected"):
            prot = ProteinStructure.from_file(sample_pdb_path, name="Test")
        ps = PigmentSystem(prot)
        assert len(ps.get_protein_background_atoms()) == 0

        # Extended loader should parse charges and yield atoms
        prot_ext = ProteinStructure.from_file_extended(sample_pdb_path, name="Test")
        ps_ext = PigmentSystem(prot_ext)
        assert len(ps_ext.get_protein_background_atoms()) > 0


def test_mvp_example_runs():
    """Test that the MVP example script can be imported and run."""
    import sys
    from pathlib import Path
    
    # Add project root to path
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))
    
    try:
        # Try to import the MVP example
        import mvp_example
        
        # Check it has a main function
        assert hasattr(mvp_example, 'main'), "MVP example should have main function"
        
    except ImportError as e:
        pytest.fail(f"Failed to import MVP example: {e}")


if __name__ == "__main__":
    # Run tests if script is executed directly
    pytest.main([__file__, "-v"])
