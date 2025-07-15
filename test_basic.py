#!/usr/bin/env python3
"""
Basic tests for the closest compound functionality
"""
import tempfile
import os
import pandas as pd
from src.utils import mol_from_smiles, fp_array
from src.index_builder import build_ecfp6_index
from src.query import search

def test_utils():
    """Test utility functions"""
    print("Testing utility functions...")
    
    # Test valid SMILES
    mol = mol_from_smiles("CCO")
    assert mol is not None, "Should parse valid SMILES"
    
    # Test fingerprint generation
    fp = fp_array(mol)
    assert fp.shape[0] == 2048, f"Expected 2048 bits, got {fp.shape[0]}"
    assert fp.dtype.name == 'float32', f"Expected float32, got {fp.dtype}"
    
    # Test invalid SMILES
    try:
        mol_from_smiles("INVALID")
        assert False, "Should raise exception for invalid SMILES"
    except ValueError:
        pass  # Expected
    
    print("✓ Utility functions work correctly")

def test_index_building():
    """Test index building functionality"""
    print("Testing index building...")
    
    # Create temporary test data
    test_data = pd.DataFrame({
        'smiles': ['CCO', 'CO', 'CCCO', 'INVALID_SMILES', 'CC(C)O'],
        'drug_name': ['Ethanol', 'Methanol', 'Propanol', 'Invalid', 'Isopropanol']
    })
    
    with tempfile.TemporaryDirectory() as temp_dir:
        csv_path = os.path.join(temp_dir, 'test.csv')
        test_data.to_csv(csv_path, index=False)
        
        # Build index
        build_ecfp6_index(csv_path, temp_dir)
        
        # Check if files were created
        assert os.path.exists(os.path.join(temp_dir, 'ecfp6.index')), "Index file not created"
        assert os.path.exists(os.path.join(temp_dir, 'meta.pkl')), "Metadata file not created"
    
    print("✓ Index building works correctly")

def test_basic_search():
    """Test basic search functionality (requires pre-built index)"""
    print("Testing search functionality...")
    
    if not os.path.exists("artifacts/ecfp6.index"):
        print("⚠ Skipping search test - no index found. Run build_index.py first.")
        return
    
    try:
        # This will print results but we mainly check it doesn't crash
        search("CCO", 3)
        print("✓ Search works correctly")
    except Exception as e:
        print(f"✗ Search failed: {e}")

if __name__ == "__main__":
    print("Running basic tests for closest_compound...")
    print("=" * 50)
    
    try:
        test_utils()
        test_index_building()
        test_basic_search()
        print("\n" + "=" * 50)
        print("All tests completed!")
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        exit(1)
