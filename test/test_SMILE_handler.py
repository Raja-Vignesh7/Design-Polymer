"""
Unit tests for SMILE_handler.py module.

Tests the following functions:
- is_valid_smiles(): SMILES validation
- draw_smile_2D(): 2D visualization and canonical SMILES
- get_smile_features(): Feature extraction from SMILES
"""

import unittest
import sys
import os
import numpy as np

# Add src/utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'utils'))

from SMILE_handler import is_valid_smiles, draw_smile_2D, get_smile_features


class TestIsValidSmiles(unittest.TestCase):
    """Test cases for is_valid_smiles() function."""
    
    def test_valid_simple_smiles(self):
        """Test validation of simple valid SMILES strings."""
        valid_smiles = [
            "C",  # Methane
            "CCO",  # Ethanol
            "c1ccccc1",  # Benzene (aromatic)
            "C1=CC=CC=C1",  # Benzene (aliphatic notation)
            "CC(C)C",  # Isobutane
        ]
        for smile in valid_smiles:
            with self.subTest(smile=smile):
                self.assertTrue(is_valid_smiles(smile))
    
    def test_valid_complex_smiles(self):
        """Test validation of complex valid SMILES strings."""
        valid_complex = [
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "c1ccccc1CC(=O)O",  # Phenylacetic acid
        ]
        for smile in valid_complex:
            with self.subTest(smile=smile):
                self.assertTrue(is_valid_smiles(smile))
    
    def test_invalid_smiles(self):
        """Test rejection of invalid SMILES strings."""
        invalid_smiles = [
            "invalid_smiles",  # Gibberish
            "C=C=C=C=",  # Valency error
            "C1CCC",  # Unclosed ring
            "Hello123!",  # Completely invalid
        ]
        for smile in invalid_smiles:
            with self.subTest(smile=smile):
                self.assertFalse(is_valid_smiles(smile))
    
    def test_empty_and_none(self):
        """Test handling of empty strings and None."""
        self.assertFalse(is_valid_smiles(""))
        self.assertFalse(is_valid_smiles(None))
    
    def test_non_string_input(self):
        """Test handling of non-string inputs."""
        self.assertFalse(is_valid_smiles(123))
        self.assertFalse(is_valid_smiles([]))
        self.assertFalse(is_valid_smiles({}))


class TestDrawSmile2D(unittest.TestCase):
    """Test cases for draw_smile_2D() function."""
    
    def test_valid_smile_returns_image_and_canonical(self):
        """Test that valid SMILES returns image and canonical SMILES."""
        smile = "CCO"
        img, canonical = draw_smile_2D(smile)
        
        self.assertIsNotNone(img)
        self.assertIsNotNone(canonical)
        self.assertEqual(canonical, "CCO")
    
    def test_canonical_smiles_output(self):
        """Test that different representations produce same canonical SMILES."""
        # Both representations should produce canonical "c1ccccc1"
        smile1 = "c1ccccc1"
        smile2 = "C1=CC=CC=C1"
        
        _, canonical1 = draw_smile_2D(smile1)
        _, canonical2 = draw_smile_2D(smile2)
        
        # Both should produce the same canonical form
        self.assertEqual(canonical1, canonical2)
    
    def test_invalid_smile_returns_none(self):
        """Test that invalid SMILES returns (None, None)."""
        img, canonical = draw_smile_2D("invalid")
        
        self.assertIsNone(img)
        self.assertIsNone(canonical)
    
    def test_image_object_type(self):
        """Test that returned image is PIL Image object."""
        smile = "CCO"
        img, _ = draw_smile_2D(smile)
        
        # Check if it has PIL Image attributes
        self.assertTrue(hasattr(img, 'size'))
        self.assertTrue(hasattr(img, 'format'))


class TestGetSmileFeatures(unittest.TestCase):
    """Test cases for get_smile_features() function."""
    
    def test_valid_smile_returns_array(self):
        """Test that valid SMILES returns numpy array."""
        smile = "CCO"
        features = get_smile_features(smile)
        
        self.assertIsInstance(features, np.ndarray)
        self.assertEqual(features.ndim, 1)
        self.assertEqual(features.shape[0], 1613)
    
    def test_invalid_smile_returns_zero_array(self):
        """Test that invalid SMILES returns zero array."""
        smile = "invalid_smiles"
        features = get_smile_features(smile)
        
        self.assertIsInstance(features, np.ndarray)
        self.assertTrue(np.all(features == 0))
        self.assertEqual(features.shape[0], 1613)
    
    def test_custom_nbits_parameter(self):
        """Test custom nBits parameter."""
        smile = "CCO"
        custom_nbits = 100
        features = get_smile_features(smile, nBits=custom_nbits)
        
        # Should return array with custom size on error
        # (valid SMILES should still return full features)
        self.assertIsInstance(features, np.ndarray)
    
    def test_features_are_numeric(self):
        """Test that all extracted features are numeric."""
        smile = "CCO"
        features = get_smile_features(smile)
        
        # All values should be numeric
        self.assertTrue(np.issubdtype(features.dtype, np.number))
    
    def test_different_smiles_different_features(self):
        """Test that different molecules produce different feature vectors."""
        smile1 = "CCO"  # Ethanol
        smile2 = "C1=CC=CC=C1"  # Benzene
        
        features1 = get_smile_features(smile1)
        features2 = get_smile_features(smile2)
        
        # Features should be different
        self.assertFalse(np.array_equal(features1, features2))
    
    def test_same_smile_produces_consistent_features(self):
        """Test that same SMILES produces consistent features."""
        smile = "CCO"
        
        features1 = get_smile_features(smile)
        features2 = get_smile_features(smile)
        
        # Features should be identical
        np.testing.assert_array_equal(features1, features2)


class TestSmileHandlerIntegration(unittest.TestCase):
    """Integration tests for SMILE_handler functions."""
    
    def test_workflow_valid_smile(self):
        """Test complete workflow with valid SMILES."""
        smile = "CCO"
        
        # Test validation
        self.assertTrue(is_valid_smiles(smile))
        
        # Test visualization
        img, canonical = draw_smile_2D(smile)
        self.assertIsNotNone(img)
        self.assertIsNotNone(canonical)
        
        # Test features
        features = get_smile_features(smile)
        self.assertEqual(features.shape[0], 1613)
    
    def test_workflow_invalid_smile(self):
        """Test complete workflow with invalid SMILES."""
        smile = "invalid_smiles"
        
        # Test validation
        self.assertFalse(is_valid_smiles(smile))
        
        # Test visualization
        img, canonical = draw_smile_2D(smile)
        self.assertIsNone(img)
        self.assertIsNone(canonical)
        
        # Test features
        features = get_smile_features(smile)
        self.assertTrue(np.all(features == 0))


if __name__ == '__main__':
    unittest.main()
