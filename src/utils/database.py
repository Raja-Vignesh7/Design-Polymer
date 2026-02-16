"""
MySQL Database handler for the Design Polymer application.
Manages connections, data retrieval, and CRUD operations.
"""

import mysql.connector
from mysql.connector import Error
import json
from typing import List, Dict, Optional, Tuple
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DatabaseHandler:
    """Handles all MySQL database operations for the Design Polymer application."""
    
    def __init__(self, host: str, user: str, password: str, database: str, port: int = 3306):
        """
        Initialize database connection parameters.
        
        Args:
            host (str): MySQL server hostname
            user (str): MySQL username
            password (str): MySQL password
            database (str): Database name
            port (int): MySQL port (default: 3306)
        """
        self.host = host
        self.user = user
        self.password = password
        self.database = database
        self.port = port
        self.connection = None
    
    def connect(self) -> bool:
        """
        Establish connection to MySQL database.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            self.connection = mysql.connector.connect(
                host=self.host,
                user=self.user,
                password=self.password,
                database=self.database,
                port=self.port
            )
            if self.connection.is_connected():
                logger.info(f"Successfully connected to database: {self.database}")
                return True
        except Error as e:
            logger.error(f"Error connecting to MySQL: {e}")
            return False
    
    def disconnect(self) -> bool:
        """
        Close database connection.
        
        Returns:
            bool: True if disconnection successful
        """
        if self.connection and self.connection.is_connected():
            self.connection.close()
            logger.info("Disconnected from MySQL database")
            return True
        return False
    
    def execute_query(self, query: str, params: Optional[Tuple] = None) -> bool:
        """
        Execute INSERT, UPDATE, or DELETE query.
        
        Args:
            query (str): SQL query to execute
            params (Optional[Tuple]): Query parameters for prepared statements
            
        Returns:
            bool: True if execution successful
        """
        try:
            cursor = self.connection.cursor()
            if params:
                cursor.execute(query, params)
            else:
                cursor.execute(query)
            self.connection.commit()
            logger.info(f"Query executed successfully. Rows affected: {cursor.rowcount}")
            cursor.close()
            return True
        except Error as e:
            logger.error(f"Error executing query: {e}")
            self.connection.rollback()
            return False
    
    def fetch_query(self, query: str, params: Optional[Tuple] = None) -> Optional[List[Dict]]:
        """
        Execute SELECT query and return results.
        
        Args:
            query (str): SQL SELECT query
            params (Optional[Tuple]): Query parameters for prepared statements
            
        Returns:
            Optional[List[Dict]]: List of dictionaries with results, None on error
        """
        try:
            cursor = self.connection.cursor(dictionary=True)
            if params:
                cursor.execute(query, params)
            else:
                cursor.execute(query)
            results = cursor.fetchall()
            cursor.close()
            return results
        except Error as e:
            logger.error(f"Error fetching data: {e}")
            return None
    
    def save_polymer_data(self, smiles: str, features: Dict, tg: float, ffv: float, 
                         density: float, tc: float, rg: float) -> bool:
        """
        Save polymer data to the database.
        
        Args:
            smiles (str): SMILES string
            features (Dict): Molecular features
            tg (float): Glass transition temperature
            ffv (float): Fractional free volume
            density (float): Polymer density
            tc (float): Thermal conductivity
            rg (float): Radius of gyration
            
        Returns:
            bool: True if data saved successfully
        """
        import numpy as np
        
        query = """
        INSERT INTO Polymer_data (SMILES, features, Tg, FFV, Density, Tc, Rg)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        
        # Sanitize features to handle inf/NaN values that aren't JSON-serializable
        def sanitize_value(val):
            """Convert inf and NaN to None (null in JSON)"""
            if isinstance(val, dict):
                return {k: sanitize_value(v) for k, v in val.items()}
            elif isinstance(val, (list, tuple)):
                return [sanitize_value(v) for v in val]
            elif isinstance(val, (float, np.floating)):
                if np.isinf(val) or np.isnan(val):
                    return float(0)  # Replace inf/NaN with 0 for JSON serialization
                return float(val)
            return val
        
        # Convert features to JSON-serializable format
        if isinstance(features, dict):
            features_clean = sanitize_value(features)
            features_json = json.dumps(features_clean)
        else:
            # Handle numpy arrays and lists
            if isinstance(features, np.ndarray):
                features_clean = sanitize_value(features.tolist())
            elif isinstance(features, (list, tuple)):
                features_clean = sanitize_value(list(features))
            else:
                features_clean = sanitize_value(features)
            features_json = json.dumps(features_clean)
        
        params = (smiles, features_json, tg, ffv, density, tc, rg)
        return self.execute_query(query, params)
    
    def get_polymer_data(self, polymer_id: Optional[int] = None) -> Optional[List[Dict]]:
        """
        Retrieve polymer data from database.
        
        Args:
            polymer_id (Optional[int]): Specific polymer ID to retrieve. If None, returns all.
            
        Returns:
            Optional[List[Dict]]: List of polymer records
        """
        if polymer_id:
            query = "SELECT * FROM Polymer_data WHERE id = %s"
            return self.fetch_query(query, (polymer_id,))
        else:
            query = "SELECT * FROM Polymer_data"
            return self.fetch_query(query)
    
    def update_polymer_data(self, polymer_id: int, **kwargs) -> bool:
        """
        Update specific polymer record.
        
        Args:
            polymer_id (int): ID of polymer to update
            **kwargs: Fields to update (e.g., Tg=100, FFV=0.25)
            
        Returns:
            bool: True if update successful
        """
        if not kwargs:
            logger.warning("No fields to update")
            return False
        
        allowed_fields = {'SMILES', 'features', 'Tg', 'FFV', 'Density', 'Tc', 'Rg'}
        
        # Filter to only allowed fields
        update_fields = {k: v for k, v in kwargs.items() if k in allowed_fields}
        
        if not update_fields:
            logger.warning("No valid fields to update")
            return False
        
        set_clause = ", ".join([f"{field} = %s" for field in update_fields.keys()])
        query = f"UPDATE Polymer_data SET {set_clause} WHERE id = %s"
        params = tuple(update_fields.values()) + (polymer_id,)
        
        return self.execute_query(query, params)
    
    def delete_polymer_data(self, polymer_id: int) -> bool:
        """
        Delete polymer record from database.
        
        Args:
            polymer_id (int): ID of polymer to delete
            
        Returns:
            bool: True if deletion successful
        """
        query = "DELETE FROM Polymer_data WHERE id = %s"
        return self.execute_query(query, (polymer_id,))
    
    def search_polymer_by_smiles(self, smiles: str) -> Optional[List[Dict]]:
        """
        Search for polymer by SMILES string.
        
        Args:
            smiles (str): SMILES string to search for
            
        Returns:
            Optional[List[Dict]]: Matching records
        """
        query = "SELECT * FROM Polymer_data WHERE SMILES = %s"
        return self.fetch_query(query, (smiles,))
    
    def save_model_metrics(self, model_type: str, version: int, mae: float, 
                          mse: float, rmse: float, r2: float, evs: float) -> bool:
        """
        Save model performance metrics.
        
        Args:
            model_type (str): Type of model (Tg, FFV, Tc, Density, Rg)
            version (int): Model version
            mae (float): Mean Absolute Error
            mse (float): Mean Squared Error
            rmse (float): Root Mean Squared Error
            r2 (float): R-squared score
            evs (float): Explained Variance Score
            
        Returns:
            bool: True if metrics saved successfully
        """
        table_map = {
            'Tg': 'Tg_models',
            'FFV': 'FFV_models',
            'Tc': 'Tc_models',
            'Density': 'Density_models',
            'Rg': 'Rg_models'
        }
        
        if model_type not in table_map:
            logger.error(f"Invalid model type: {model_type}")
            return False
        
        table = table_map[model_type]
        query = f"""
        INSERT INTO {table} (version, MAE, MSE, RMSE, R2, EVS)
        VALUES (%s, %s, %s, %s, %s, %s)
        """
        params = (version, mae, mse, rmse, r2, evs)
        return self.execute_query(query, params)
    
    def get_model_metrics(self, model_type: str, version: Optional[int] = None) -> Optional[List[Dict]]:
        """
        Retrieve model metrics.
        
        Args:
            model_type (str): Type of model (Tg, FFV, Tc, Density, Rg)
            version (Optional[int]): Specific version to retrieve. If None, returns all.
            
        Returns:
            Optional[List[Dict]]: Model metrics
        """
        table_map = {
            'Tg': 'Tg_models',
            'FFV': 'FFV_models',
            'Tc': 'Tc_models',
            'Density': 'Density_models',
            'Rg': 'Rg_models'
        }
        
        if model_type not in table_map:
            logger.error(f"Invalid model type: {model_type}")
            return None
        
        table = table_map[model_type]
        
        if version:
            query = f"SELECT * FROM {table} WHERE version = %s"
            return self.fetch_query(query, (version,))
        else:
            query = f"SELECT * FROM {table} ORDER BY version DESC"
            return self.fetch_query(query)
    
    def get_statistics(self) -> Optional[Dict]:
        """
        Get database statistics.
        
        Returns:
            Optional[Dict]: Statistics about polymers and models
        """
        try:
            cursor = self.connection.cursor(dictionary=True)
            
            stats = {}
            
            # Count polymers
            cursor.execute("SELECT COUNT(*) as count FROM Polymer_data")
            stats['total_polymers'] = cursor.fetchone()['count']
            
            # Average properties
            cursor.execute("""
            SELECT AVG(Tg) as avg_tg, AVG(FFV) as avg_ffv, AVG(Density) as avg_density,
                   AVG(Tc) as avg_tc, AVG(Rg) as avg_rg
            FROM Polymer_data
            """)
            averages = cursor.fetchone()
            stats['averages'] = averages
            
            # Model counts
            models = ['Tg_models', 'FFV_models', 'Tc_models', 'Density_models', 'Rg_models']
            for model in models:
                cursor.execute(f"SELECT COUNT(*) as count FROM {model}")
                stats[f'{model}_count'] = cursor.fetchone()['count']
            
            cursor.close()
            return stats
        except Error as e:
            logger.error(f"Error getting statistics: {e}")
            return None


# # Example usage functions
# def example_usage():
#     """
#     Example of how to use the DatabaseHandler class.
#     """
#     # Initialize database handler
#     db = DatabaseHandler(
#         host='localhost',
#         user='root',
#         password='root',
#         database='polymerdb'
#     )
    
#     # Connect to database
#     if db.connect():
#         # Save polymer data
#         features = {'molecular_weight': 100, 'logp': 2.5}
#         db.save_polymer_data(
#             smiles='C1=CC=CC=C1',
#             features=features,
#             tg=120.5,
#             ffv=0.28,
#             density=1.2,
#             tc=0.15,
#             rg=5.3
#         )
        
#         # Retrieve all polymers
#         polymers = db.get_polymer_data()
#         print(polymers)
        
#         # Update polymer data
#         db.update_polymer_data(1, Tg=125.0, FFV=0.30)
        
#         # Save model metrics
#         db.save_model_metrics(
#             model_type='Tg',
#             version=1,
#             mae=5.2,
#             mse=27.0,
#             rmse=5.2,
#             r2=0.95,
#             evs=0.94
#         )
        
#         # Get statistics
#         stats = db.get_statistics()
#         print(stats)
        
#         # Disconnect
#         db.disconnect()


# if __name__ == "__main__":
#     example_usage()
