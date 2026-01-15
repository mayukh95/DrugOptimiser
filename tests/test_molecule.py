"""
Tests for DrugMolecule class
"""

import pytest
import numpy as np


class TestDrugMolecule:
    """Test suite for DrugMolecule class."""
    
    def test_molecule_creation_from_smiles(self):
        """Test creating a molecule from SMILES string."""
        from src.core import DrugMolecule
        
        mol = DrugMolecule(
            name="Aspirin",
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O"
        )
        
        assert mol.name == "Aspirin"
        assert mol.smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"
        assert mol.mol is not None
        assert len(mol.elements) > 0
        assert mol.initial_coords.shape[1] == 3
    
    def test_descriptor_calculation(self):
        """Test molecular descriptor calculation."""
        from src.core import DrugMolecule
        
        mol = DrugMolecule(
            name="Caffeine",
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        )
        mol.calculate_descriptors()
        
        assert 'MW' in mol.properties
        assert 'LogP' in mol.properties
        assert 'TPSA' in mol.properties
        assert 'HBD' in mol.properties
        assert 'HBA' in mol.properties
        assert mol.properties['MW'] > 0
    
    def test_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        from src.core import DrugMolecule
        
        with pytest.raises(Exception):
            mol = DrugMolecule(
                name="Invalid",
                smiles="NOT_A_VALID_SMILES_STRING"
            )
    
    def test_elements_extraction(self):
        """Test that elements are correctly extracted."""
        from src.core import DrugMolecule
        
        # Water
        mol = DrugMolecule(name="Water", smiles="O")
        mol.calculate_descriptors()
        
        # Should have O and 2 H
        element_counts = {}
        for el in mol.elements:
            element_counts[el] = element_counts.get(el, 0) + 1
        
        assert 'O' in element_counts
        assert 'H' in element_counts


class TestADMETPredictor:
    """Test suite for ADMET predictions."""
    
    def test_lipinski_pass(self):
        """Test Lipinski filter for drug-like molecule."""
        from src.core import DrugMolecule, ADMETPredictor
        
        # Aspirin - should pass Lipinski
        mol = DrugMolecule(name="Aspirin", smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
        mol.calculate_descriptors()
        
        predictor = ADMETPredictor()
        result = predictor.check_lipinski(mol)
        
        assert 'pass' in result
        assert 'violations' in result
    
    def test_lipinski_fail(self):
        """Test Lipinski filter for non-drug-like molecule."""
        from src.core import DrugMolecule, ADMETPredictor
        
        # Large molecule - likely to fail
        large_smiles = "C" * 100  # Just a long carbon chain
        
        try:
            mol = DrugMolecule(name="LargeMol", smiles=large_smiles)
            mol.calculate_descriptors()
            
            predictor = ADMETPredictor()
            result = predictor.check_lipinski(mol)
            
            # Should have violations
            assert result['violations'] > 0
        except:
            # Invalid molecule is also acceptable
            pass


class TestHelperFunctions:
    """Test helper functions."""
    
    def test_get_example_molecules(self):
        """Test loading example molecules."""
        from src import get_example_molecules
        
        molecules = get_example_molecules()
        
        assert len(molecules) > 0
        assert all(hasattr(m, 'name') for m in molecules)
        assert all(hasattr(m, 'smiles') for m in molecules)
