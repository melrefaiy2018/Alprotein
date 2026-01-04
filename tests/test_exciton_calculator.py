#!/usr/bin/env python3
"""
Test ExcitonCalculator integration.

This script tests the ExcitonCalculator class to ensure it correctly
integrates with the AlProtein workflow.
"""

import numpy as np
import pandas as pd
import tempfile
from pathlib import Path
import sys

# Add Alprotein to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from Alprotein.calculators.exciton_calculator import ExcitonCalculator


def create_test_hamiltonian():
    """Create a simple test Hamiltonian."""
    # Create a small test system with 4 pigments
    pigment_names = ['A_CLA_501', 'A_CLA_502', 'A_CLA_503', 'A_CLA_504']

    # Simple Hamiltonian with diagonal energies and some couplings
    H = np.array([
        [15000, 50, 10, 5],
        [50, 14900, 30, 8],
        [10, 30, 14950, 40],
        [5, 8, 40, 15100]
    ], dtype=float)

    H_df = pd.DataFrame(H, index=pigment_names, columns=pigment_names)
    return H_df


def create_test_domains():
    """Create test domain structure."""
    # Two domains: {0: [0, 1], 1: [2, 3]}
    return {
        0: [0, 1],  # First two pigments
        1: [2, 3]   # Last two pigments
    }


def create_mock_pigment_system():
    """Create a mock PigmentSystem for testing."""
    class MockPigment:
        def __init__(self, name):
            self.name = name
            self.dipole = np.array([0.0, 0.0, 1.0])  # Simple z-direction dipole

    class MockPigmentSystem:
        def __init__(self):
            self.dict_pigments = {
                'A_CLA_501': MockPigment('A_CLA_501'),
                'A_CLA_502': MockPigment('A_CLA_502'),
                'A_CLA_503': MockPigment('A_CLA_503'),
                'A_CLA_504': MockPigment('A_CLA_504')
            }

    return MockPigmentSystem()


def test_initialization():
    """Test ExcitonCalculator initialization."""
    print("\n[TEST 1] Testing initialization...")

    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=100, temperature=300.0)

    assert calc.disorder_sigma == 50.0, "Disorder sigma not set correctly"
    assert calc.n_ensemble == 100, "Ensemble size not set correctly"
    assert calc.temperature == 300.0, "Temperature not set correctly"
    assert calc.distributions is None, "Distributions should be None before calculation"
    assert calc.site_labels is None, "Site labels should be None before calculation"

    print("    ✓ Initialization test passed")
    return calc


def test_domain_conversion():
    """Test domain format conversion from indices to names."""
    print("\n[TEST 2] Testing domain conversion...")

    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=10)

    domains = {0: [0, 1], 1: [2, 3]}
    pigment_names = ['A_CLA_501', 'A_CLA_502', 'A_CLA_503', 'A_CLA_504']

    converted = calc._convert_domains_to_names(domains, pigment_names)

    # Check structure
    assert isinstance(converted, list), "Output should be a list"
    assert len(converted) == 2, "Should have 2 domains"

    # Check first domain
    assert converted[0] == ['A_CLA_501', 'A_CLA_502'], "First domain incorrect"

    # Check second domain
    assert converted[1] == ['A_CLA_503', 'A_CLA_504'], "Second domain incorrect"

    print("    ✓ Domain conversion test passed")
    print(f"      Input:  {domains}")
    print(f"      Output: {converted}")


def test_calculate_distributions():
    """Test exciton distribution calculation."""
    print("\n[TEST 3] Testing calculate_distributions...")

    # Create test data
    hamiltonian = create_test_hamiltonian()
    domains = create_test_domains()
    pigment_system = create_mock_pigment_system()

    # Create calculator with small ensemble for speed
    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=10)

    # Calculate distributions
    distributions, labels = calc.calculate_distributions(
        hamiltonian=hamiltonian,
        domains=domains,
        pigment_system=pigment_system
    )

    # Verify output structure
    assert isinstance(distributions, dict), "Distributions should be a dict"
    assert isinstance(labels, list), "Labels should be a list"
    assert len(labels) == 4, "Should have 4 pigment labels"

    # Check that all pigments have distributions
    for pigment in labels:
        assert pigment in distributions, f"Missing distribution for {pigment}"
        energies, probs = distributions[pigment]
        assert isinstance(energies, list), "Energies should be a list"
        assert isinstance(probs, list), "Probabilities should be a list"
        assert len(energies) == len(probs), "Energies and probs should have same length"
        assert len(energies) > 0, "Should have at least one energy value"

    # Check that results are stored
    assert calc.distributions is not None, "Distributions should be stored"
    assert calc.site_labels is not None, "Site labels should be stored"

    print("    ✓ Calculate distributions test passed")
    print(f"      Calculated distributions for {len(labels)} pigments")
    print(f"      Example (A_CLA_501): {len(distributions['A_CLA_501'][0])} data points")


def test_error_handling():
    """Test error handling when plotting before calculation."""
    print("\n[TEST 4] Testing error handling...")

    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=10)

    # Try to plot before calculating
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            calc.plot_distributions(output_path=tmpdir)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "No distributions calculated" in str(e)
        print("    ✓ Correctly raised ValueError for plot_distributions")

    # Try combined plot before calculating
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            calc.plot_combined_with_absorption(
                wavelengths_abs=np.linspace(600, 720, 100),
                absorption=np.random.rand(100),
                output_path=tmpdir
            )
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "No distributions calculated" in str(e)
        print("    ✓ Correctly raised ValueError for plot_combined_with_absorption")

    # Try export before calculating
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            calc.export_distributions(str(Path(tmpdir) / "test.csv"))
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "No distributions calculated" in str(e)
        print("    ✓ Correctly raised ValueError for export_distributions")


def test_csv_export():
    """Test CSV export functionality."""
    print("\n[TEST 5] Testing CSV export...")

    # Create test data and calculate
    hamiltonian = create_test_hamiltonian()
    domains = create_test_domains()
    pigment_system = create_mock_pigment_system()

    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=10)
    calc.calculate_distributions(hamiltonian, domains, pigment_system)

    # Export to temporary file
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = Path(tmpdir) / "exciton_distributions.csv"
        calc.export_distributions(str(csv_path))

        # Verify file was created
        assert csv_path.exists(), "CSV file should be created"

        # Read and verify contents
        df = pd.read_csv(csv_path)

        # Check columns
        expected_columns = {'pigment', 'energy_cm-1', 'wavelength_nm', 'probability'}
        assert set(df.columns) == expected_columns, f"CSV should have columns {expected_columns}"

        # Check that we have data for all pigments
        pigments_in_csv = set(df['pigment'].unique())
        expected_pigments = {'A_CLA_501', 'A_CLA_502', 'A_CLA_503', 'A_CLA_504'}
        assert pigments_in_csv == expected_pigments, "Should have all pigments in CSV"

        # Check that energies and wavelengths are consistent
        for _, row in df.head(10).iterrows():
            if row['energy_cm-1'] > 0 and not np.isnan(row['wavelength_nm']):
                expected_wl = 1e7 / row['energy_cm-1']
                assert abs(row['wavelength_nm'] - expected_wl) < 0.1, "Wavelength conversion incorrect"

        print("    ✓ CSV export test passed")
        print(f"      Exported {len(df)} rows")
        print(f"      Pigments: {len(pigments_in_csv)}")


def test_repr():
    """Test string representation."""
    print("\n[TEST 6] Testing __repr__...")

    calc = ExcitonCalculator(disorder_sigma=55.3, n_ensemble=300, temperature=300.0)
    repr_str = repr(calc)

    assert "ExcitonCalculator" in repr_str
    assert "55.3" in repr_str
    assert "300" in repr_str

    print(f"    ✓ __repr__ test passed")
    print(f"      {repr_str}")


def test_integration_workflow():
    """Test full integration workflow."""
    print("\n[TEST 7] Testing full integration workflow...")

    # Create test data
    hamiltonian = create_test_hamiltonian()
    domains = create_test_domains()
    pigment_system = create_mock_pigment_system()

    # Initialize calculator (use small ensemble for speed)
    calc = ExcitonCalculator(disorder_sigma=50.0, n_ensemble=20)

    print("    - Calculating distributions...")
    distributions, labels = calc.calculate_distributions(
        hamiltonian=hamiltonian,
        domains=domains,
        pigment_system=pigment_system
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Export CSV
        print("    - Exporting to CSV...")
        csv_path = tmppath / "distributions.csv"
        calc.export_distributions(str(csv_path))
        assert csv_path.exists()

        # Note: We skip plotting tests here since they require matplotlib
        # and display, which may not be available in CI/CD environments.
        # The plotting functions delegate to proven existing code anyway.

    print("    ✓ Full integration workflow test passed")


def run_all_tests():
    """Run all tests."""
    print("="*80)
    print("EXCITON CALCULATOR TESTS")
    print("="*80)

    try:
        test_initialization()
        test_domain_conversion()
        test_calculate_distributions()
        test_error_handling()
        test_csv_export()
        test_repr()
        test_integration_workflow()

        print("\n" + "="*80)
        print("✓ ALL TESTS PASSED!")
        print("="*80)
        return True

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
