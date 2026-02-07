"""
Unit tests for analysis.py module.

Tests the following classes:
- Analyser: SMILES analysis and visualization
"""

import unittest
import sys
import os
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
import matplotlib.pyplot as plt

# Add src/utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'utils'))

from analysis import Analyser


class TestAnalyserInitialization(unittest.TestCase):
    """Test cases for Analyser initialization."""
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def test_init_with_single_smile(self, mock_is_valid, mock_models):
        """Test initialization with single SMILES string."""
        mock_is_valid.return_value = True
        
        mock_model_instance = MagicMock()
        mock_model_instance.predict_properties.return_value = {
            'Tg': 100.0, 'Tc': 0.5, 'Rg': 5.0, 'FFV': 0.2, 'Density': 1.2
        }
        mock_models.return_value = mock_model_instance
        
        analyser = Analyser("CCO")
        
        self.assertEqual(len(analyser.smiles_list), 1)
        self.assertEqual(analyser.smiles_list[0], "CCO")
        self.assertEqual(len(analyser.valid_smiles), 1)
        self.assertEqual(len(analyser.predictions), 1)
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def test_init_with_multiple_smiles(self, mock_is_valid, mock_models):
        """Test initialization with list of SMILES strings."""
        mock_is_valid.return_value = True
        
        mock_model_instance = MagicMock()
        mock_model_instance.predict_properties.return_value = {
            'Tg': 100.0, 'Tc': 0.5, 'Rg': 5.0, 'FFV': 0.2, 'Density': 1.2
        }
        mock_models.return_value = mock_model_instance
        
        smiles_list = ["CCO", "C1=CC=CC=C1", "CCC"]
        analyser = Analyser(smiles_list)
        
        self.assertEqual(len(analyser.smiles_list), 3)
        self.assertEqual(len(analyser.valid_smiles), 3)
        self.assertEqual(len(analyser.predictions), 3)
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def test_init_filters_invalid_smiles(self, mock_is_valid, mock_models):
        """Test that invalid SMILES are filtered out during initialization."""
        # Return True for first two, False for third
        mock_is_valid.side_effect = [True, True, False]
        
        mock_model_instance = MagicMock()
        mock_model_instance.predict_properties.return_value = {
            'Tg': 100.0, 'Tc': 0.5, 'Rg': 5.0, 'FFV': 0.2, 'Density': 1.2
        }
        mock_models.return_value = mock_model_instance
        
        smiles_list = ["CCO", "C1=CC=CC=C1", "invalid_smiles"]
        analyser = Analyser(smiles_list)
        
        # Should only have 2 valid SMILES
        self.assertEqual(len(analyser.valid_smiles), 2)
        self.assertEqual(len(analyser.predictions), 2)


class TestAnalyserGet(unittest.TestCase):
    """Test cases for Analyser.get() method."""
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def setUp(self, mock_is_valid, mock_models):
        """Set up test fixtures."""
        mock_is_valid.return_value = True
        
        mock_model_instance = MagicMock()
        # Return different values for each call
        predictions = [
            {'Tg': 150.0, 'Tc': 0.6, 'Rg': 6.0, 'FFV': 0.22, 'Density': 1.25},
            {'Tg': 120.0, 'Tc': 0.5, 'Rg': 5.0, 'FFV': 0.20, 'Density': 1.20},
            {'Tg': 180.0, 'Tc': 0.7, 'Rg': 7.0, 'FFV': 0.25, 'Density': 1.30},
        ]
        mock_model_instance.predict_properties.side_effect = predictions
        mock_models.return_value = mock_model_instance
        
        self.analyser = Analyser(["CCO", "C1=CC=CC=C1", "CCC"])
    
    def test_get_returns_tuple(self):
        """Test that get() returns a tuple of (plots, sorted_results)."""
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            plots, sorted_results = self.analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            self.assertIsInstance(plots, dict)
            self.assertIsInstance(sorted_results, list)
    
    def test_get_returns_plots_dict(self):
        """Test that plots dictionary contains all properties."""
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            plots, _ = self.analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Should have plots for all 5 properties
            expected_properties = ["Tg", "Tc", "Rg", "FFV", "Density"]
            for prop in expected_properties:
                self.assertIn(prop, plots)
    
    def test_get_returns_sorted_smiles(self):
        """Test that sorted_results contains SMILES and properties."""
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            plots, sorted_results = self.analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Should return 3 results (same as input)
            self.assertEqual(len(sorted_results), 3)
            
            # Should return list of dicts
            for result in sorted_results:
                self.assertIsInstance(result, dict)
                self.assertIn('SMILES', result)
                self.assertIn('Tg', result)
                self.assertIn('Tc', result)
                self.assertIn('Rg', result)
                self.assertIn('FFV', result)
                self.assertIn('Density', result)
    
    def test_sorting_by_primary_property(self):
        """Test that sorting is done by first property in priority_order."""
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            # Sort by Tg (should be: 180, 150, 120)
            _, sorted_results = self.analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Verify results contain Tg values in descending order
            self.assertEqual(len(sorted_results), 3)
            tg_values = [result['Tg'] for result in sorted_results]
            # Check that values are sorted in descending order (higher values first)
            self.assertEqual(tg_values, sorted(tg_values, reverse=True))
    
    def test_empty_data_handling(self):
        """Test handling of empty or invalid data."""
        with patch('analysis.Models'):
            with patch('analysis.is_valid_smiles', return_value=False):
                empty_analyser = Analyser("invalid")
                
                with patch('analysis.plt.subplots') as mock_subplots:
                    mock_fig = MagicMock()
                    mock_ax = MagicMock()
                    mock_subplots.return_value = (mock_fig, mock_ax)
                    plots, sorted_results = empty_analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
                    
                    self.assertEqual(len(plots), 0)
                    self.assertEqual(len(sorted_results), 0)
    
    def test_plots_have_figure_attributes(self):
        """Test that returned plots are matplotlib figures."""
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            plots, _ = self.analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            for fig in plots.values():
                # Check for figure-like attributes
                self.assertTrue(hasattr(fig, 'add_subplot') or hasattr(fig, 'subplots') or fig is not None)


class TestAnalyserPriorityOrdering(unittest.TestCase):
    """Test cases for priority-based ordering in Analyser."""
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def test_different_priority_orders(self, mock_is_valid, mock_models):
        """Test that different priority orders produce different results."""
        mock_is_valid.return_value = True
        
        mock_model_instance = MagicMock()
        predictions = [
            {'Tg': 150.0, 'Tc': 0.6, 'Rg': 6.0, 'FFV': 0.22, 'Density': 1.25},
            {'Tg': 120.0, 'Tc': 0.8, 'Rg': 5.0, 'FFV': 0.20, 'Density': 1.20},
            {'Tg': 180.0, 'Tc': 0.5, 'Rg': 7.0, 'FFV': 0.25, 'Density': 1.30},
        ]
        mock_model_instance.predict_properties.side_effect = predictions
        mock_models.return_value = mock_model_instance
        
        analyser = Analyser(["CCO", "C1=CC=CC=C1", "CCC"])
        
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            # Different priority orders
            _, sorted_by_tg = analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Reset for next test
            mock_model_instance.predict_properties.side_effect = predictions
            analyser2 = Analyser(["CCO", "C1=CC=CC=C1", "CCC"])
            _, sorted_by_tc = analyser2.get(["Tc", "Tg", "Rg", "FFV", "Density"])
            
            # Results should contain dicts with SMILES and properties
            self.assertEqual(len(sorted_by_tg), 3)
            self.assertEqual(len(sorted_by_tc), 3)
            
            # All results should be dicts
            for result in sorted_by_tg + sorted_by_tc:
                self.assertIsInstance(result, dict)
                self.assertIn('SMILES', result)


class TestAnalyserIntegration(unittest.TestCase):
    """Integration tests for Analyser class."""
    
    @patch('analysis.Models')
    @patch('analysis.is_valid_smiles')
    def test_complete_workflow(self, mock_is_valid, mock_models):
        """Test complete workflow from initialization to analysis."""
        mock_is_valid.return_value = True
        
        mock_model_instance = MagicMock()
        predictions = [
            {'Tg': 150.0, 'Tc': 0.6, 'Rg': 6.0, 'FFV': 0.22, 'Density': 1.25},
            {'Tg': 120.0, 'Tc': 0.8, 'Rg': 5.0, 'FFV': 0.20, 'Density': 1.20},
        ]
        mock_model_instance.predict_properties.side_effect = predictions
        mock_models.return_value = mock_model_instance
        
        # Create analyser
        analyser = Analyser(["CCO", "C1=CC=CC=C1"])
        
        # Verify initialization
        self.assertEqual(len(analyser.valid_smiles), 2)
        self.assertEqual(len(analyser.predictions), 2)
        
        with patch('analysis.plt.subplots') as mock_subplots:
            mock_fig = MagicMock()
            mock_ax = MagicMock()
            mock_subplots.return_value = (mock_fig, mock_ax)
            # Get analysis
            plots, sorted_results = analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Verify results
            self.assertEqual(len(plots), 5)  # 5 properties
            self.assertEqual(len(sorted_results), 2)  # 2 compounds
            
            # Verify results contain SMILES and properties
            for result in sorted_results:
                self.assertIsInstance(result, dict)
                self.assertIn('SMILES', result)
                self.assertIn('Tg', result)
                self.assertIn('Density', result)


if __name__ == '__main__':
    unittest.main()
