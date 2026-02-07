from dotenv import load_dotenv
import warnings
import os
import numpy as np
from SMILE_handler import get_smile_features
import warnings
import joblib
warnings.filterwarnings('ignore')

class db_config_info:
    """
    Database configuration manager for handling .env file operations.
    
    This class manages the database configuration stored in a .env file,
    allowing retrieval and setting of database connection parameters.
    """
    
    def __init__(self):
        """
        Initialize the configuration manager.
        
        Loads environment variables and checks for the existence of .env file.
        Creates the file if it doesn't exist.
        """
        load_dotenv()
        self.FILE_PATH = 'src\\utils\\.env'
        
        # Create .env file if it doesn't exist
        if not os.path.exists(self.FILE_PATH):
            with open(self.FILE_PATH, 'w') as f:
                pass
        
        # Check if .env file contains any configuration
        with open(self.FILE_PATH, 'r') as f:
            lines = f.readlines()
            print(lines)
            self.contains_default = len(lines) > 0
    def get_config_info(self):
        """
        Retrieve database configuration from .env file.
        
        Returns:
            dict: Dictionary containing database configuration with keys:
                  DB_HOST, DB_USER, DB_PASSWORD, DB_NAME, DB_PORT
                  
        Raises:
            Exception: If .env file is not found or failed to load configuration
        """
        if self.contains_default:
            # Load environment variables from .env file
            if load_dotenv(dotenv_path=self.FILE_PATH):
                config_info = {
                    "DB_HOST": os.getenv("DB_HOST"),
                    "DB_USER": os.getenv("DB_USER"),
                    "DB_PASSWORD": os.getenv("DB_PASSWORD"),
                    "DB_NAME": os.getenv("DB_NAME"),
                    "DB_PORT": os.getenv("DB_PORT")
                }
                return config_info
            else:
                raise Exception("Failed to load .env file.")
        else:
            raise Exception("No Database configuration found. Please set up the .env file.")
    
    def set_config_info(self, host, user, password, db_name, port):
        """
        Write database configuration to .env file.
        
        Args:
            host (str): Database host address
            user (str): Database username
            password (str): Database password
            db_name (str): Database name
            port (int or str): Database port number
        """
        with open(self.FILE_PATH, 'w') as f:
            f.write(f"DB_HOST={host}\n")
            f.write(f"DB_USER={user}\n")
            f.write(f"DB_PASSWORD={password}\n")
            f.write(f"DB_NAME={db_name}\n")
            f.write(f"DB_PORT={port}\n")
        self.contains_default = True
class Models:
    """
    Machine Learning Models Manager for polymer property prediction.
    
    This class loads pre-trained machine learning models for predicting various
    polymer properties (Tg, Tc, Rg, FFV, Density) from SMILES strings.
    
    Models are loaded from joblib-serialized files for efficient computation.
    """
    
    def __init__(self):
        """
        Initialize and load all pre-trained ML models.
        
        Loads five separate joblib model files for predicting different polymer properties:
        - Tg_model: Glass transition temperature prediction
        - Tc_model: Thermal conductivity prediction
        - Rg_model: Radius of gyration prediction
        - FFV_model: Fractional free volume prediction
        - Density_model: Polymer density prediction
        
        Raises:
            FileNotFoundError: If any model file is not found at the specified path
        """
        # Define model directory path relative to current script location
        model_dir = os.path.join('..', 'models')
        
        # Load pre-trained Tg (Glass Transition Temperature) model
        self.tg_model = joblib.load( 'src\\models\\Tg_model.joblib')
        
        # Load pre-trained Tc (Thermal Conductivity) model
        self.Tc_model = joblib.load( 'src\\models\\Tc_model.joblib')
        
        # Load pre-trained Rg (Radius of Gyration) model
        self.Rg_model = joblib.load( 'src\\models\\Rg_model.joblib')
        
        # Load pre-trained FFV (Fractional Free Volume) model
        self.FFV_model = joblib.load( 'src\\models\\FFV_model.joblib')
        
        # Load pre-trained Density model
        self.Density_model = joblib.load( 'src\\models\\Density_model.joblib')
        
    def predict_properties(self, smile):
        """
        Predict all polymer properties from SMILES string(s).
        
        This method extracts molecular features from the input SMILES string(s) and
        uses the pre-trained models to predict multiple polymer properties in one call.
        
        Args:
            smile (str or list): SMILES string(s) representing polymer fragment(s).
                                Can be a single string or list of strings.
            
        Returns:
            dict or list: If input is a single SMILES string, returns a dictionary 
                         containing predicted properties. If input is a list, returns 
                         a list of dictionaries for each SMILES.
                         
                         Dictionary keys:
                         - "Tg" (float): Predicted glass transition temperature
                         - "Tc" (float): Predicted thermal conductivity
                         - "Rg" (float): Predicted radius of gyration
                         - "FFV" (float): Predicted fractional free volume
                         - "Density" (float): Predicted polymer density
                  
        Example:
            >>> models = Models()
            >>> # Single SMILES
            >>> predictions = models.predict_properties("CCO")
            >>> print(f"Predicted Tg: {predictions['Tg']}")
            >>> 
            >>> # Multiple SMILES
            >>> predictions = models.predict_properties(["CCO", "C1=CC=CC=C1"])
            >>> for pred in predictions:
            ...     print(f"Predicted Tg: {pred['Tg']}")
            
        Notes:
            - Feature extraction uses Mordred descriptors (1800 features)
            - All models expect the same feature vector format
            - Returns single float values (not arrays)
            - Batch processing for lists is more efficient than individual calls
        """
        # Check if input is a list of SMILES strings
        if isinstance(smile, list):
            # Extract features for all SMILES strings
            features_list = [get_smile_features(smi) for smi in smile]
            features_array = np.array(features_list)
            
            # Predict all properties for all SMILES at once (batch prediction)
            tg_preds = self.tg_model.predict(features_array)
            Tc_preds = self.Tc_model.predict(features_array)
            Rg_preds = self.Rg_model.predict(features_array)
            FFV_preds = self.FFV_model.predict(features_array)
            Density_preds = self.Density_model.predict(features_array)
            
            # Format results as list of dictionaries
            results = []
            for i in range(len(smile)):
                results.append({
                    "Tg": float(tg_preds[i]),
                    "Tc": float(Tc_preds[i]),
                    "Rg": float(Rg_preds[i]),
                    "FFV": float(FFV_preds[i]),
                    "Density": float(Density_preds[i])
                })
            return results
        
        # Handle single SMILES string
        # Extract molecular features from SMILES using Mordred descriptors
        features = get_smile_features(smile)
        
        # Predict Tg using the trained model
        # Takes feature array and returns single prediction value
        tg_pred = self.tg_model.predict([features])[0]
        
        # Predict Tc (Thermal Conductivity) using the trained model
        Tc_pred = self.Tc_model.predict([features])[0]
        
        # Predict Rg (Radius of Gyration) using the trained model
        Rg_pred = self.Rg_model.predict([features])[0]
        
        # Predict FFV (Fractional Free Volume) using the trained model
        FFV_pred = self.FFV_model.predict([features])[0]
        
        # Predict Density using the trained model
        Density_pred = self.Density_model.predict([features])[0]
        
        # Return all predictions as a dictionary for easy access
        return {
            "Tg": tg_pred,
            "Tc": Tc_pred,
            "Rg": Rg_pred,
            "FFV": FFV_pred,
            "Density": Density_pred
        }
    
# db_configer = db_config_info()
# print(db_configer.get_config_info())
# ip = "*CC(*)c1ccccc1C(=O)OCCCCCC"  # Benzene
# model = Models()
# print(model.predict_properties(ip))