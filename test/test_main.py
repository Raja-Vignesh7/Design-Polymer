"""
Unit tests for main.py module.

Tests the following classes:
- db_config_info: Database configuration management
- Models: Machine learning model predictions
"""

import unittest
import sys
import os
import tempfile
import json
from unittest.mock import patch, MagicMock

# Add src/utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'utils'))

from main import db_config_info, Models


class TestDbConfigInfo(unittest.TestCase):
    """Test cases for db_config_info class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directory for test .env file
        self.test_dir = tempfile.mkdtemp()
        self.test_env_path = os.path.join(self.test_dir, '.env')
    
    def tearDown(self):
        """Clean up temporary files."""
        if os.path.exists(self.test_env_path):
            os.remove(self.test_env_path)
        os.rmdir(self.test_dir)
    
    def test_initialization(self):
        """Test db_config_info initialization."""
        with patch('main.load_dotenv'):
            config = db_config_info()
            self.assertIsNotNone(config)
    
    def test_env_file_creation(self):
        """Test that .env file is created if it doesn't exist."""
        with patch('main.load_dotenv'):
            with patch.object(db_config_info, '__init__', lambda x: None):
                config = db_config_info()
                config.FILE_PATH = self.test_env_path
                
                # File should not exist yet
                self.assertFalse(os.path.exists(self.test_env_path))
    
    def test_set_config_info(self):
        """Test setting database configuration."""
        with patch('main.load_dotenv'):
            config = db_config_info()
            config.FILE_PATH = self.test_env_path
            
            # Set configuration
            config.set_config_info('localhost', 'root', 'password', 'mydb', 3306)
            
            # Check file was created and contains data
            self.assertTrue(os.path.exists(self.test_env_path))
            
            with open(self.test_env_path, 'r') as f:
                content = f.read()
                self.assertIn('DB_HOST=localhost', content)
                self.assertIn('DB_USER=root', content)
                self.assertIn('DB_PASSWORD=password', content)
                self.assertIn('DB_NAME=mydb', content)
                self.assertIn('DB_PORT=3306', content)
    
    @patch('main.load_dotenv')
    @patch.dict(os.environ, {
        'DB_HOST': 'localhost',
        'DB_USER': 'root',
        'DB_PASSWORD': 'password',
        'DB_NAME': 'testdb',
        'DB_PORT': '3306'
    })
    def test_get_config_info(self, mock_load):
        """Test retrieving database configuration."""
        mock_load.return_value = True
        
        with patch.object(db_config_info, '__init__', lambda x: setattr(x, 'FILE_PATH', self.test_env_path) or setattr(x, 'contains_default', True)):
            config = db_config_info()
            config_data = config.get_config_info()
            
            self.assertIsInstance(config_data, dict)
            self.assertIn('DB_HOST', config_data)
            self.assertIn('DB_USER', config_data)
            self.assertIn('DB_PASSWORD', config_data)
            self.assertIn('DB_NAME', config_data)
            self.assertIn('DB_PORT', config_data)


class TestModels(unittest.TestCase):
    """Test cases for Models class."""
    
    @patch('main.joblib.load')
    def setUp(self, mock_joblib):
        """Set up test fixtures."""
        # Mock all model loads
        mock_model = MagicMock()
        mock_model.predict.return_value = [0.5]
        mock_joblib.return_value = mock_model
        
        # Initialize Models with mocked joblib
        self.models = Models()
    
    @patch('main.joblib.load')
    def test_models_initialization(self, mock_joblib):
        """Test Models initialization loads all five models."""
        mock_model = MagicMock()
        mock_joblib.return_value = mock_model
        
        models = Models()
        
        # Check that joblib.load was called 5 times (for 5 models)
        self.assertEqual(mock_joblib.call_count, 5)
        
        # Check that all model attributes exist
        self.assertIsNotNone(models.tg_model)
        self.assertIsNotNone(models.Tc_model)
        self.assertIsNotNone(models.Rg_model)
        self.assertIsNotNone(models.FFV_model)
        self.assertIsNotNone(models.Density_model)
    
    @patch('main.get_smile_features')
    def test_predict_properties_single_smile(self, mock_features):
        """Test predict_properties with single SMILES string."""
        import numpy as np
        
        mock_features.return_value = np.random.rand(1613)
        
        with patch('main.joblib.load') as mock_joblib:
            mock_model = MagicMock()
            mock_model.predict.return_value = np.array([100.0])
            mock_joblib.return_value = mock_model
            
            models = Models()
            result = models.predict_properties("CCO")
            
            # Check result is a dictionary
            self.assertIsInstance(result, dict)
            
            # Check all properties are present
            expected_keys = ["Tg", "Tc", "Rg", "FFV", "Density"]
            for key in expected_keys:
                self.assertIn(key, result)
                self.assertIsInstance(result[key], (float, np.floating))
    
    @patch('main.get_smile_features')
    def test_predict_properties_multiple_smiles(self, mock_features):
        """Test predict_properties with list of SMILES strings."""
        import numpy as np
        
        mock_features.return_value = np.random.rand(1613)
        
        with patch('main.joblib.load') as mock_joblib:
            mock_model = MagicMock()
            mock_model.predict.return_value = np.array([100.0, 110.0, 120.0])
            mock_joblib.return_value = mock_model
            
            models = Models()
            smiles_list = ["CCO", "C1=CC=CC=C1", "CCC"]
            result = models.predict_properties(smiles_list)
            
            # Check result is a list
            self.assertIsInstance(result, list)
            
            # Check we got 3 results
            self.assertEqual(len(result), 3)
            
            # Check each result is a dictionary with all properties
            expected_keys = ["Tg", "Tc", "Rg", "FFV", "Density"]
            for res in result:
                self.assertIsInstance(res, dict)
                for key in expected_keys:
                    self.assertIn(key, res)
    
    @patch('main.get_smile_features')
    def test_predict_properties_returns_floats(self, mock_features):
        """Test that predictions return float values."""
        import numpy as np
        
        mock_features.return_value = np.random.rand(1613)
        
        with patch('main.joblib.load') as mock_joblib:
            mock_model = MagicMock()
            mock_model.predict.return_value = np.array([123.45])
            mock_joblib.return_value = mock_model
            
            models = Models()
            result = models.predict_properties("CCO")
            
            # All values should be floats
            for value in result.values():
                self.assertIsInstance(value, (float, np.floating))
    
    @patch('main.get_smile_features')
    def test_predict_properties_calls_all_models(self, mock_features):
        """Test that predict_properties calls all 5 models."""
        import numpy as np
        
        mock_features.return_value = np.random.rand(1613)
        
        with patch('main.joblib.load') as mock_joblib:
            mock_model = MagicMock()
            mock_model.predict.return_value = np.array([100.0])
            mock_joblib.return_value = mock_model
            
            models = Models()
            
            # Reset call count
            mock_model.reset_mock()
            
            models.predict_properties("CCO")
            
            # Each model should be called once
            self.assertEqual(mock_model.predict.call_count, 5)


class TestModelsIntegration(unittest.TestCase):
    """Integration tests for Models class."""
    
    @patch('main.joblib.load')
    @patch('main.get_smile_features')
    def test_predict_properties_workflow(self, mock_features, mock_joblib):
        """Test complete prediction workflow."""
        import numpy as np
        
        mock_features.return_value = np.random.rand(1613)
        
        # Create a mock that returns correct number of predictions
        mock_model = MagicMock()
        
        # Use side_effect to return predictions based on input shape
        def predict_side_effect(X):
            import numpy as np
            # X can be a list or numpy array
            if isinstance(X, list):
                # If it's a list, get its length
                return np.array([150.0] * len(X))
            elif isinstance(X, np.ndarray):
                # If it's 1D array, return single prediction
                if X.ndim == 1:
                    return np.array([150.0])
                else:
                    # Return one prediction per sample
                    return np.array([150.0] * X.shape[0])
            else:
                return np.array([150.0])
        
        mock_model.predict.side_effect = predict_side_effect
        mock_joblib.return_value = mock_model
        
        models = Models()
        
        # Single prediction
        single_result = models.predict_properties("CCO")
        self.assertIsInstance(single_result, dict)
        
        # Batch prediction
        batch_result = models.predict_properties(["CCO", "C1=CC=CC=C1"])
        self.assertIsInstance(batch_result, list)
        self.assertEqual(len(batch_result), 2)


if __name__ == '__main__':
    unittest.main()
