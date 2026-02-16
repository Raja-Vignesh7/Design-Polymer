"""
SMILE_handler.py - SMILES String Processing and Molecular Feature Extraction

This module provides utilities for validating, visualizing, and extracting molecular features
from SMILES (Simplified Molecular Input Line Entry System) strings using RDKit and Mordred.

Key functionalities:
- Validate SMILES format and chemical correctness
- Generate 2D molecular visualizations
- Extract molecular descriptors and features using Mordred
- Convert SMILES to feature vectors for machine learning

Dependencies:
- rdkit: RDKit chemistry toolkit for molecular operations
- mordred: Molecular descriptor calculator
- numpy: Numerical computing
- matplotlib: Visualization (optional)
"""

# Import required libraries
from rdkit import Chem  # RDKit core module for molecular operations
from rdkit.Chem import AllChem  # RDKit module for advanced chemical operations
from rdkit import rdBase  # RDKit base module
from mordred import Calculator, descriptors  # Mordred for molecular descriptor calculation
from rdkit.Chem import Draw  # RDKit module for drawing molecules
import numpy as np  # Numerical computing library
import matplotlib.pyplot as plt  # Plotting library


def is_valid_smiles(ip):
    """
    Validates if a given input string is a valid SMILES format.
    
    This function checks both the syntactic validity of the SMILES string and ensures
    that the molecular structure is chemically valid (correct valency, ring closure, etc.).
    
    Args:
        ip (str): Input string to validate as SMILES format
        
    Returns:
        bool: True if the input is a valid SMILES string, False otherwise
        
    Example:
        >>> is_valid_smiles("C1=CC=CC=C1")  # Benzene
        True
        >>> is_valid_smiles("invalid_smiles")
        False
        
    Notes:
        - Returns False if input is None, empty, or not a string
        - MolFromSmiles with sanitize=True ensures chemical correctness
    """
    # Validate input type and non-empty check
    if not ip or not isinstance(ip, str):
        return False
        
    # MolFromSmiles returns a Molecule object if valid, None if invalid
    # sanitize=True ensures the chemistry (valency, atomic bonding, etc.) is also correct
    mol = Chem.MolFromSmiles(ip, sanitize=True)
    
    return mol is not None


def draw_smile_2D(smile):
    """
    Generates a 2D visualization of a molecular structure from SMILES and returns canonical SMILES.
    
    This function creates a visual representation of the molecule and standardizes the SMILES
    representation to a canonical form (consistent ordering of atoms and bonds).
    
    Args:
        smile (str): SMILES string representing a molecular structure
        
    Returns:
        tuple: (img, canonical_smiles) where:
            - img: PIL Image object containing the 2D molecular structure visualization
            - canonical_smiles (str): Canonical form of the SMILES string
            Returns (None, None) if the input SMILES is invalid
            
    Example:
        >>> img, canon_smiles = draw_smile_2D("CCO")  # Ethanol
        >>> # img can be displayed in Streamlit or saved to file
        >>> print(canon_smiles)  # 'CCO'
        
    Notes:
        - The function validates SMILES before processing
        - Canonical SMILES is useful for database lookups and comparisons
        - Image output can be directly used with Streamlit's st.image()
    """
    # Validate SMILES format before processing
    if not is_valid_smiles(smile):
        return None, None
    
    # Convert SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smile)
    
    # Generate 2D image representation of the molecule
    img = Draw.MolToImage(mol)
    
    # Convert molecule to canonical SMILES (standardized representation)
    # Canonical form ensures consistent string representation regardless of input format
    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
    
    return img, canon_smiles


def get_smile_features(smile,as_dict=False, nBits=1613):
    """
    Extracts molecular descriptors from a SMILES string and returns them as a feature vector.
    
    This function calculates ~1613 molecular descriptors using Mordred, which includes
    topological, geometric, and physicochemical properties. Invalid values are replaced with 0.
    
    Args:
        smile (str): SMILES string representing a molecular structure
        nBits (int): Expected number of features/descriptors (default: 1613)
                    Used as array size for invalid SMILES
        
    Returns:
        numpy.ndarray: 1D numpy array of molecular features/descriptors
                      - Shape: (1613,) by default
                      - Returns array of zeros if SMILES is invalid or feature extraction fails
                      - Non-numeric or invalid descriptor values are converted to 0.0
            
    Example:
        >>> features = get_smile_features("C1=CC=CC=C1")  # Benzene
        >>> print(features.shape)  # (1613,)
        >>> print(features[:5])  # First 5 features
        
    Notes:
        - Uses Mordred calculator configured for 2D descriptors (ignore_3D=True)
        - Robust error handling: returns zero vector on any failure
        - Suitable for machine learning applications (neural networks, regression models)
        - Processing time depends on molecule size and complexity
    """
    try:
        # Validate SMILES format first
        if not is_valid_smiles(smile):
            # Return zero vector for invalid SMILES
            return np.zeros(nBits)
        
        # Convert SMILES string to RDKit molecule object
        mol = Chem.MolFromSmiles(smile)
        
        # Initialize Mordred descriptor calculator
        # ignore_3D=True means we only use 2D descriptors (faster computation)
        calc = Calculator(descriptors, ignore_3D=True)
        
        # Calculate all molecular descriptors for the molecule
        mordred_desc = calc(mol)
        
        # Convert descriptor object to dictionary format
        # Keys are descriptor names, values are descriptor values
        mordred_df = mordred_desc.asdict()
        if as_dict:
            for key in mordred_df:
                try:
                    mordred_df[key] = float(mordred_df[key])
                except (ValueError, TypeError):
                    mordred_df[key] = 0.0
            return mordred_df
        # Initialize list to store feature values
        arr = []
        
        # Iterate through all descriptor values
        for val in mordred_df.values():
            try:
                # Attempt to convert descriptor value to float
                if int(val):  # Check if value can be converted to int (non-zero)
                    arr.append(float(val))
                else:
                    # Zero values are preserved as 0
                    arr.append(0)
            except (ValueError, TypeError):
                # Handle non-numeric descriptor values by replacing with 0
                arr.append(0)
        
        # Convert list to numpy array for efficient computation
        result = np.array(arr)
        
        # Return the actual features (no padding/truncation)
        # Mordred typically returns ~1613 features which matches trained models
        return result
        
    except Exception as e:
        # Catch-all error handling: log exception and return zero vector
        print(f"Error extracting features from SMILES '{smile}': {e}")
        return np.zeros(nBits)


# ============================================================================
# EXAMPLE USAGE AND TEST CASES
# ============================================================================
# The following examples demonstrate how to use the SMILE_handler functions
# Uncomment to test or run the examples below
# ============================================================================

# Example 1: Validate SMILES strings
# test_cases = {
#     "Benzene": "c1ccccc1",                              # Valid aromatic ring
#     "Ethanol": "CCO",                                   # Valid simple alcohol
#     "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",        # Valid complex molecule
#     "Invalid Carbon": "C=C=C=C=",                       # Invalid (Valency issue)
#     "Broken Ring": "C1CCC",                             # Invalid (Unclosed ring)
#     "Gibberish": "Hello123!"                            # Invalid (Not SMILES)
# }
#
# print(f"{'Name':<15} | {'SMILES':<30} | {'Valid?':<7}")
# print("-" * 55)
# for name, smiles in test_cases.items():
#     valid = is_valid_smiles(smiles)
#     print(f"{name:<15} | {smiles:<30} | {str(valid):<7}")


# Example 2: Visualize and get canonical SMILES
# smile_input = "C1=CC=CC=C1"  # Benzene
# 
# if is_valid_smiles(smile_input):
#     img, canonical_smile = draw_smile_2D(smile_input)
#     if img:
#         # Display in Streamlit
#         # st.image(img, caption=f"Canonical SMILES: {canonical_smile}")
#         
#         # Or display with matplotlib
#         plt.imshow(img)
#         plt.axis('off')
#         plt.title(f"Canonical SMILES: {canonical_smile}")
#         plt.show()
# else:
#     print("Invalid SMILES input")


# Example 3: Extract molecular features
# smile_input = "CCO"  # Ethanol
# features = get_smile_features(smile_input)
# 
# if features is not None:
#     print(f"Number of features extracted: {len(features)}")
#     print(f"Feature vector shape: {features.shape}")
#     print(f"First 10 features: {features[:10]}")
#     print(f"Non-zero features: {np.count_nonzero(features)}")
# else:
#     print("Failed to extract features")