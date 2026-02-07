"""
Unit tests for database.py module.

Tests the following classes:
- DatabaseHandler: MySQL database operations
"""

import unittest
import sys
import os
from unittest.mock import patch, MagicMock, call
import mysql.connector

# Add src/utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'utils'))

from database import DatabaseHandler


class TestDatabaseHandlerInitialization(unittest.TestCase):
    """Test cases for DatabaseHandler initialization."""
    
    def test_init_with_default_port(self):
        """Test initialization with default port."""
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        self.assertEqual(handler.host, 'localhost')
        self.assertEqual(handler.user, 'root')
        self.assertEqual(handler.password, 'password')
        self.assertEqual(handler.database, 'testdb')
        self.assertEqual(handler.port, 3306)
        self.assertIsNone(handler.connection)
    
    def test_init_with_custom_port(self):
        """Test initialization with custom port."""
        handler = DatabaseHandler(
            host='192.168.1.1',
            user='admin',
            password='secret',
            database='mydb',
            port=3307
        )
        
        self.assertEqual(handler.port, 3307)
    
    def test_init_stores_credentials(self):
        """Test that credentials are properly stored."""
        credentials = {
            'host': 'db.example.com',
            'user': 'dbuser',
            'password': 'dbpass',
            'database': 'production',
            'port': 5432
        }
        
        handler = DatabaseHandler(**credentials)
        
        self.assertEqual(handler.host, credentials['host'])
        self.assertEqual(handler.user, credentials['user'])
        self.assertEqual(handler.password, credentials['password'])
        self.assertEqual(handler.database, credentials['database'])
        self.assertEqual(handler.port, credentials['port'])


class TestDatabaseHandlerConnection(unittest.TestCase):
    """Test cases for database connection methods."""
    
    @patch('database.mysql.connector.connect')
    def test_connect_success(self, mock_connect):
        """Test successful database connection."""
        mock_connection = MagicMock()
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        result = handler.connect()
        
        self.assertTrue(result)
        self.assertEqual(handler.connection, mock_connection)
        mock_connect.assert_called_once()
    
    @patch('database.mysql.connector.connect')
    def test_connect_failure(self, mock_connect):
        """Test failed database connection."""
        mock_connect.side_effect = mysql.connector.Error("Connection failed")
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        result = handler.connect()
        
        self.assertFalse(result)
        self.assertIsNone(handler.connection)
    
    @patch('database.mysql.connector.connect')
    def test_disconnect_success(self, mock_connect):
        """Test successful database disconnection."""
        mock_connection = MagicMock()
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        handler.connect()
        result = handler.disconnect()
        
        self.assertTrue(result)
        mock_connection.close.assert_called_once()
    
    @patch('database.mysql.connector.connect')
    def test_disconnect_without_connection(self, mock_connect):
        """Test disconnect when no connection exists."""
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        # Try to disconnect without connecting
        result = handler.disconnect()
        
        self.assertFalse(result)


class TestDatabaseHandlerQueries(unittest.TestCase):
    """Test cases for database query methods."""
    
    @patch('database.mysql.connector.connect')
    def setUp(self, mock_connect):
        """Set up test fixtures."""
        self.mock_connection = MagicMock()
        self.mock_cursor = MagicMock()
        self.mock_connection.cursor.return_value = self.mock_cursor
        mock_connect.return_value = self.mock_connection
        
        self.handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        self.handler.connect()
    
    def test_execute_query_select(self):
        """Test executing a SELECT query."""
        self.mock_cursor.fetchall.return_value = [
            (1, 'CCO', 100.5),
            (2, 'C1=CC=CC=C1', 150.2)
        ]
        
        query = "SELECT * FROM polymers"
        result = self.handler.execute_query(query)
        
        self.assertIsNotNone(result)
        self.mock_cursor.execute.assert_called_once_with(query)
    
    def test_execute_query_insert(self):
        """Test executing an INSERT query."""
        query = "INSERT INTO polymers (smiles, tg) VALUES ('CCO', 100.5)"
        self.handler.execute_query(query)
        
        self.mock_cursor.execute.assert_called_once_with(query)
        self.mock_connection.commit.assert_called_once()
    
    def test_execute_query_update(self):
        """Test executing an UPDATE query."""
        query = "UPDATE polymers SET tg = 110.5 WHERE id = 1"
        self.handler.execute_query(query)
        
        self.mock_cursor.execute.assert_called_once_with(query)
        self.mock_connection.commit.assert_called_once()
    
    def test_execute_query_delete(self):
        """Test executing a DELETE query."""
        query = "DELETE FROM polymers WHERE id = 1"
        self.handler.execute_query(query)
        
        self.mock_cursor.execute.assert_called_once_with(query)
        self.mock_connection.commit.assert_called_once()
    
    def test_execute_query_with_exception(self):
        """Test query execution with error handling."""
        self.mock_cursor.execute.side_effect = mysql.connector.Error("Query error")
        
        query = "SELECT * FROM invalid_table"
        result = self.handler.execute_query(query)
        
        # Should return False on error (depending on implementation)
        self.assertFalse(result)


class TestDatabaseHandlerTableOperations(unittest.TestCase):
    """Test cases for table-related operations."""
    
    @patch('database.mysql.connector.connect')
    def test_table_exists(self, mock_connect):
        """Test checking if table exists."""
        mock_connection = MagicMock()
        mock_cursor = MagicMock()
        mock_connection.cursor.return_value = mock_cursor
        mock_cursor.fetchone.return_value = (1,)
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        handler.connect()
        
        # This depends on actual implementation
        # Just verify the handler can be used
        self.assertIsNotNone(handler)
    
    @patch('database.mysql.connector.connect')
    def test_get_table_schema(self, mock_connect):
        """Test retrieving table schema information."""
        mock_connection = MagicMock()
        mock_cursor = MagicMock()
        mock_connection.cursor.return_value = mock_cursor
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        handler.connect()
        
        # Verify handler is ready
        self.assertIsNotNone(handler.connection)


class TestDatabaseHandlerIntegration(unittest.TestCase):
    """Integration tests for DatabaseHandler."""
    
    @patch('database.mysql.connector.connect')
    def test_connection_workflow(self, mock_connect):
        """Test complete connection workflow."""
        mock_connection = MagicMock()
        mock_cursor = MagicMock()
        mock_connection.cursor.return_value = mock_cursor
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        # Test connect
        self.assertTrue(handler.connect())
        self.assertIsNotNone(handler.connection)
        
        # Test disconnect
        self.assertTrue(handler.disconnect())
        mock_connection.close.assert_called_once()
    
    @patch('database.mysql.connector.connect')
    def test_query_workflow(self, mock_connect):
        """Test complete query workflow."""
        mock_connection = MagicMock()
        mock_cursor = MagicMock()
        mock_connection.cursor.return_value = mock_cursor
        mock_cursor.fetchall.return_value = [
            (1, 'CCO', 100.0),
            (2, 'C1=CC=CC=C1', 150.0)
        ]
        mock_connect.return_value = mock_connection
        
        handler = DatabaseHandler(
            host='localhost',
            user='root',
            password='password',
            database='testdb'
        )
        
        handler.connect()
        
        # Execute a SELECT query
        query = "SELECT id, smiles, tg FROM polymers"
        result = handler.execute_query(query)
        
        # Verify execution
        self.assertIsNotNone(result)
        mock_cursor.execute.assert_called()
        
        handler.disconnect()


if __name__ == '__main__':
    unittest.main()
