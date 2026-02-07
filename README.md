# Design Polymer

Design Polymer is a tool that helps in creating and designing polymers with predictions for these 5 desired properties:

1. **Tg** - Glass transition temperature (°C)
2. **FFV** - Fractional free volume
3. **Tc** - Thermal conductivity (W/m·K)
4. **Density** - Polymer density (g/cm³)
5. **Rg** - Radius of gyration

## Installation

Clone the repository and install the required dependencies:

```bash
git clone https://github.com/Raja-Vignesh7/Design-Polymer.git
cd "Design Polymer"
pip install -r requirements.txt
```

## Usage

### Running the Streamlit Application

Start the interactive web application:

```bash
cd src
streamlit run app.py
```

The application will open in your default browser at `http://localhost:8501`.

### Using the Python Modules

You can also use the modules directly in your Python code:

```python
from src.utils.SMILE_handler import process_smile
from src.utils.analysis import predict_properties
```

## Project Structure

```
Design Polymer/
├── src/                          # Main application source code
│   ├── app.py                   # Streamlit web application
│   ├── models/                  # Pre-trained ML models
│   │   ├── Density_model.joblib
│   │   ├── FFV_model.joblib
│   │   ├── Rg_model.joblib
│   │   ├── Tc_model.joblib
│   │   └── Tg_model.joblib
│   └── utils/                   # Utility modules
│       ├── analysis.py          # Property analysis functions
│       ├── database.py          # Database operations
│       ├── main.py              # Main logic
│       └── SMILE_handler.py     # SMILES notation processing
├── test/                         # Unit tests
│   ├── test_analysis.py
│   ├── test_app.py
│   ├── test_database.py
│   ├── test_main.py
│   └── test_SMILE_handler.py
├── SQL queries/                  # Database SQL scripts
│   └── create tables.sql
├── README.md                     # This file
└── requirements.txt              # Python dependencies
```

## Features

- **SMILES Input Support** - Process polymer structures using SMILES notation
- **Multi-Property Prediction** - Predict 5 key polymer properties simultaneously
- **Machine Learning Models** - Pre-trained models for accurate predictions
- **Interactive Web Interface** - User-friendly Streamlit application
- **Database Integration** - Store and retrieve polymer data
- **Comprehensive Testing** - Full test coverage for all modules

## Requirements

- Python 3.8+
- Dependencies listed in `requirements.txt`

Key packages:
- streamlit - Web application framework
- scikit-learn - Machine learning library
- pandas - Data manipulation
- numpy - Numerical computing

## Testing

Run the test suite using Python's unittest framework:

```bash
# Run all tests
python -m unittest discover -s test -p "test_*.py" -v

# Run specific test module
python -m unittest test.test_analysis -v
```

All tests should pass successfully.
