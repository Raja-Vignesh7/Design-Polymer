# Test Suite Documentation

This directory contains comprehensive unit and integration tests for the Design Polymer application.

## Test Files

### 1. **test_SMILE_handler.py**
Tests for SMILES processing and molecular feature extraction:
- `is_valid_smiles()` - SMILES string validation
- `draw_smile_2D()` - 2D molecular visualization
- `get_smile_features()` - Feature extraction using Mordred

**Test Classes:**
- `TestIsValidSmiles` - Validation logic tests
- `TestDrawSmile2D` - Visualization tests
- `TestGetSmileFeatures` - Feature extraction tests
- `TestSmileHandlerIntegration` - Integration tests

### 2. **test_main.py**
Tests for machine learning models and configuration management:
- `db_config_info` - Database configuration handling
- `Models` - ML model predictions

**Test Classes:**
- `TestDbConfigInfo` - Configuration management tests
- `TestModels` - Model prediction tests
- `TestModelsIntegration` - Model workflow tests

### 3. **test_analysis.py**
Tests for the Analyser class that generates visualizations and sorted results:
- `Analyser.__init__()` - Initialization and prediction
- `Analyser.get()` - Analysis and visualization generation

**Test Classes:**
- `TestAnalyserInitialization` - Initialization tests
- `TestAnalyserGet` - Analysis and sorting tests
- `TestAnalyserPriorityOrdering` - Priority-based sorting tests
- `TestAnalyserIntegration` - Complete workflow tests

### 4. **test_database.py**
Tests for MySQL database operations:
- `DatabaseHandler` - Database connection and CRUD operations

**Test Classes:**
- `TestDatabaseHandlerInitialization` - Connection setup tests
- `TestDatabaseHandlerConnection` - Connect/disconnect tests
- `TestDatabaseHandlerQueries` - Query execution tests
- `TestDatabaseHandlerTableOperations` - Table operation tests
- `TestDatabaseHandlerIntegration` - Complete workflow tests

### 5. **test_app.py**
Tests for the Streamlit application:
- App configuration and UI components
- Input validation
- Layout testing

**Test Classes:**
- `TestAppImports` - Module import tests
- `TestAppConfiguration` - App setup tests
- `TestAppUIComponents` - UI component tests
- `TestAppValidation` - Input validation tests
- `TestAppLayout` - Layout component tests
- `TestAppIntegration` - Complete app workflow tests

## Running Tests

### Run all tests:
```bash
python -m unittest discover -s test -p "test_*.py"
```

### Run specific test file:
```bash
python -m unittest test.test_SMILE_handler
```

### Run specific test class:
```bash
python -m unittest test.test_analysis.TestAnalyserGet
```

### Run with verbose output:
```bash
python -m unittest discover -s test -p "test_*.py" -v
```

### Run with coverage report (requires coverage):
```bash
pip install coverage
coverage run -m unittest discover -s test -p "test_*.py"
coverage report
coverage html  # Generate HTML report
```

## Test Coverage

The test suite covers:
- ✅ Function validation (valid/invalid inputs)
- ✅ Error handling and edge cases
- ✅ Data type verification
- ✅ Integration between modules
- ✅ Mocking external dependencies (ML models, database)
- ✅ Output format and structure validation

## Dependencies

Required for testing:
- `unittest` (built-in)
- `numpy` and `pandas`
- `matplotlib`
- `rdkit` and `mordred`
- `streamlit`
- `mysql-connector-python`

For coverage analysis:
- `coverage`

## Notes

1. **Mocking**: Heavy use of `unittest.mock` to mock external dependencies like database connections and ML models, allowing tests to run without requiring actual model files or database setup.

2. **SMILE_handler Tests**: Tests verify both valid and invalid SMILES strings, feature consistency, and visualization capabilities.

3. **Models Tests**: Tests use mocked joblib.load to avoid requiring actual model files (.joblib).

4. **Analyser Tests**: Tests verify sorting by priority order, plot generation, and handling of batch SMILES processing.

5. **Database Tests**: Tests use mocked MySQL connections to test query execution without a live database.

6. **Streamlit App Tests**: Tests are limited due to Streamlit's UI-based nature; focus is on module imports and basic functionality.

## Continuous Integration

These tests can be integrated into CI/CD pipelines:
- GitHub Actions
- GitLab CI
- Jenkins
- Travis CI

Example GitHub Actions workflow in `.github/workflows/tests.yml`:
```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
      - run: pip install -r requirements.txt
      - run: python -m pytest test/
```
