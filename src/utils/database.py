"""
SQLite database handler for the Design Polymer application.
Manages connections, schema creation, data retrieval, and CRUD operations.
"""

import json
import logging
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DatabaseHandler:
    """Handles SQLite database operations for the Design Polymer application."""

    def __init__(
        self,
        host: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
        database: Optional[str] = None,
        port: int = 3306,
        database_path: Optional[Union[str, Path]] = None,
    ):
        """Initialize database connection parameters."""
        self.host = host or "localhost"
        self.user = user or ""
        self.password = password or ""
        self.database = database or "polymer_db"
        self.port = port

        if database_path is not None:
            self.database_path = Path(database_path)
        else:
            database_name = str(self.database)
            self.database_path = Path(database_name if Path(database_name).suffix else f"{database_name}.db")

        self.connection: Optional[sqlite3.Connection] = None

    def _normalize_query(self, query: str) -> str:
        """Convert MySQL-style placeholders to SQLite placeholders."""
        return query.replace("%s", "?")

    def _initialize_schema(self) -> None:
        """Create required tables if they do not already exist."""
        if not self.connection:
            return

        schema = [
            """
            CREATE TABLE IF NOT EXISTS Polymer_data (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                SMILES TEXT NOT NULL,
                features TEXT NOT NULL,
                Tg REAL NOT NULL,
                FFV REAL NOT NULL,
                Density REAL NOT NULL,
                Tc REAL NOT NULL,
                Rg REAL NOT NULL
            )
            """,
            """
            CREATE TABLE IF NOT EXISTS Tg_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                version INTEGER NOT NULL,
                MAE REAL NOT NULL,
                MSE REAL NOT NULL,
                RMSE REAL NOT NULL,
                R2 REAL NOT NULL,
                EVS REAL NOT NULL
            )
            """,
            """
            CREATE TABLE IF NOT EXISTS FFV_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                version INTEGER NOT NULL,
                MAE REAL NOT NULL,
                MSE REAL NOT NULL,
                RMSE REAL NOT NULL,
                R2 REAL NOT NULL,
                EVS REAL NOT NULL
            )
            """,
            """
            CREATE TABLE IF NOT EXISTS Tc_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                version INTEGER NOT NULL,
                MAE REAL NOT NULL,
                MSE REAL NOT NULL,
                RMSE REAL NOT NULL,
                R2 REAL NOT NULL,
                EVS REAL NOT NULL
            )
            """,
            """
            CREATE TABLE IF NOT EXISTS Density_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                version INTEGER NOT NULL,
                MAE REAL NOT NULL,
                MSE REAL NOT NULL,
                RMSE REAL NOT NULL,
                R2 REAL NOT NULL,
                EVS REAL NOT NULL
            )
            """,
            """
            CREATE TABLE IF NOT EXISTS Rg_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                version INTEGER NOT NULL,
                MAE REAL NOT NULL,
                MSE REAL NOT NULL,
                RMSE REAL NOT NULL,
                R2 REAL NOT NULL,
                EVS REAL NOT NULL
            )
            """,
        ]

        for statement in schema:
            self.connection.execute(statement)
        self.connection.commit()

    def connect(self) -> bool:
        """Establish connection to the SQLite database."""
        try:
            if self.connection is None:
                self.connection = sqlite3.connect(str(self.database_path))
                self.connection.row_factory = sqlite3.Row
                self._initialize_schema()
            logger.info(f"Successfully connected to database: {self.database_path}")
            return True
        except sqlite3.Error as exc:
            logger.error(f"Error connecting to SQLite: {exc}")
            return False

    def disconnect(self) -> bool:
        """Close the database connection."""
        if self.connection:
            self.connection.close()
            self.connection = None
            logger.info("Disconnected from database")
            return True
        return False

    def execute_query(self, query: str, params: Optional[Tuple] = None) -> bool:
        """Execute an INSERT, UPDATE, or DELETE query."""
        if not self.connection:
            return False

        try:
            cursor = self.connection.cursor()
            normalized_query = self._normalize_query(query)
            if params is not None:
                cursor.execute(normalized_query, params)
            else:
                cursor.execute(normalized_query)
            self.connection.commit()
            logger.info(f"Query executed successfully. Rows affected: {cursor.rowcount}")
            cursor.close()
            return True
        except sqlite3.Error as exc:
            logger.error(f"Error executing query: {exc}")
            if self.connection:
                self.connection.rollback()
            return False

    def fetch_query(self, query: str, params: Optional[Tuple] = None) -> Optional[List[Dict]]:
        """Execute a SELECT query and return results as dictionaries."""
        if not self.connection:
            return None

        try:
            cursor = self.connection.cursor()
            normalized_query = self._normalize_query(query)
            if params is not None:
                cursor.execute(normalized_query, params)
            else:
                cursor.execute(normalized_query)
            results = cursor.fetchall()
            cursor.close()
            return [dict(row) for row in results]
        except sqlite3.Error as exc:
            logger.error(f"Error fetching data: {exc}")
            return None

    def save_polymer_data(self, smiles: str, features: Dict, tg: float, ffv: float,
                         density: float, tc: float, rg: float) -> bool:
        """Save polymer data to the database."""
        import numpy as np

        query = """
        INSERT INTO Polymer_data (SMILES, features, Tg, FFV, Density, Tc, Rg)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """

        def sanitize_value(val):
            """Convert inf and NaN to None for JSON serialization."""
            if isinstance(val, dict):
                return {k: sanitize_value(v) for k, v in val.items()}
            if isinstance(val, (list, tuple)):
                return [sanitize_value(v) for v in val]
            if isinstance(val, (float, np.floating)):
                if np.isinf(val) or np.isnan(val):
                    return float(0)
                return float(val)
            return val

        if isinstance(features, dict):
            features_clean = sanitize_value(features)
            features_json = json.dumps(features_clean)
        else:
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
        """Retrieve polymer data from the database."""
        if polymer_id is not None:
            query = "SELECT * FROM Polymer_data WHERE id = ?"
            return self.fetch_query(query, (polymer_id,))

        query = "SELECT * FROM Polymer_data"
        return self.fetch_query(query)

    def update_polymer_data(self, polymer_id: int, **kwargs) -> bool:
        """Update a specific polymer record."""
        if not kwargs:
            logger.warning("No fields to update")
            return False

        allowed_fields = {'SMILES', 'features', 'Tg', 'FFV', 'Density', 'Tc', 'Rg'}
        update_fields = {k: v for k, v in kwargs.items() if k in allowed_fields}

        if not update_fields:
            logger.warning("No valid fields to update")
            return False

        set_clause = ", ".join([f"{field} = ?" for field in update_fields.keys()])
        query = f"UPDATE Polymer_data SET {set_clause} WHERE id = ?"
        params = tuple(update_fields.values()) + (polymer_id,)

        return self.execute_query(query, params)

    def delete_polymer_data(self, polymer_id: int) -> bool:
        """Delete a polymer record from the database."""
        query = "DELETE FROM Polymer_data WHERE id = ?"
        return self.execute_query(query, (polymer_id,))

    def search_polymer_by_smiles(self, smiles: str) -> Optional[List[Dict]]:
        """Search for polymer records by SMILES string."""
        query = "SELECT * FROM Polymer_data WHERE SMILES = ?"
        return self.fetch_query(query, (smiles,))

    def save_model_metrics(self, model_type: str, version: int, mae: float,
                          mse: float, rmse: float, r2: float, evs: float) -> bool:
        """Save model performance metrics."""
        table_map = {
            'Tg': 'Tg_models',
            'FFV': 'FFV_models',
            'Tc': 'Tc_models',
            'Density': 'Density_models',
            'Rg': 'Rg_models',
        }

        if model_type not in table_map:
            logger.error(f"Invalid model type: {model_type}")
            return False

        table = table_map[model_type]
        query = f"""
        INSERT INTO {table} (version, MAE, MSE, RMSE, R2, EVS)
        VALUES (?, ?, ?, ?, ?, ?)
        """
        params = (version, mae, mse, rmse, r2, evs)
        return self.execute_query(query, params)

    def get_model_metrics(self, model_type: str, version: Optional[int] = None) -> Optional[List[Dict]]:
        """Retrieve model metrics."""
        table_map = {
            'Tg': 'Tg_models',
            'FFV': 'FFV_models',
            'Tc': 'Tc_models',
            'Density': 'Density_models',
            'Rg': 'Rg_models',
        }

        if model_type not in table_map:
            logger.error(f"Invalid model type: {model_type}")
            return None

        table = table_map[model_type]
        if version is not None:
            query = f"SELECT * FROM {table} WHERE version = ?"
            return self.fetch_query(query, (version,))

        query = f"SELECT * FROM {table} ORDER BY version DESC"
        return self.fetch_query(query)

    def get_statistics(self) -> Optional[Dict]:
        """Get database statistics."""
        if not self.connection:
            return None

        try:
            cursor = self.connection.cursor()
            stats = {}

            cursor.execute("SELECT COUNT(*) as count FROM Polymer_data")
            stats['total_polymers'] = cursor.fetchone()['count']

            cursor.execute("""
            SELECT AVG(Tg) as avg_tg, AVG(FFV) as avg_ffv, AVG(Density) as avg_density,
                   AVG(Tc) as avg_tc, AVG(Rg) as avg_rg
            FROM Polymer_data
            """)
            averages = cursor.fetchone()
            stats['averages'] = dict(averages)

            models = ['Tg_models', 'FFV_models', 'Tc_models', 'Density_models', 'Rg_models']
            for model in models:
                cursor.execute(f"SELECT COUNT(*) as count FROM {model}")
                stats[f'{model}_count'] = cursor.fetchone()['count']

            cursor.close()
            return stats
        except sqlite3.Error as exc:
            logger.error(f"Error getting statistics: {exc}")
            return None
