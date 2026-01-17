from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from mordred import Calculator, descriptors
from rdkit.Chem import Draw
import numpy as np
import matplotlib.pyplot as plt


def is_valid_smiles(ip):
    """
    Checks if a string is a valid SMILES format using RDKit.
    """
    if not ip or not isinstance(ip, str):
        return False
        
    # MolFromSmiles returns a Molecule object if valid, None if invalid
    # sanitize=True ensures the chemistry (valency, etc.) is also correct
    mol = Chem.MolFromSmiles(ip, sanitize=True)
    
    return mol is not None

def get_smile_features(smile):
    if not is_valid_smiles(smile):
        return None
    mol = Chem.MolFromSmiles(smile)
    calc = Calculator(descriptors, ignore_3D=True)
    features = calc(mol)
    return features

def draw_smile_2D(smile):
    if not is_valid_smiles(smile):
        return None, None
    mol = Chem.MolFromSmiles(smile)
    img = Draw.MolToImage(mol)
    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
    return img, canon_smiles

def get_smile_features(smile,nBits=1800):
        try:
            if not is_valid_smiles(smile):
                return np.zeros(nBits)
            mol = Chem.MolFromSmiles(smile)
            calc = Calculator(descriptors, ignore_3D=True)
            mordred_desc = calc(mol)
            
            # Convert to DataFrame
            mordred_df = mordred_desc.asdict()
            arr = []
            for val in mordred_df.values():
                try:
                    if int(val):
                        arr.append(float(val))
                    else:
                        arr.append(0)
                except:
                    arr.append(0)
            return np.array(arr)
        except Exception as e:
            print(e)
            return np.zeros(nBits)

# ip = "C1=CC=CC=C1"  # Benzene


# if not ip or not isinstance(ip, str):
#     print("Invalid input")
# elif Chem.MolFromSmiles(ip) is None:
#     print("Invalid SMILES")
# else:
#     mol = Chem.MolFromSmiles(ip)
#     print(mol)
#     canon_smiles = Chem.MolToSmiles(mol, canonical=True)

#     img = Draw.MolToImage(mol)
#     plt.imshow(img)
#     plt.axis('off')
#     plt.show()
#     print(canon_smiles)

# # --- Testing the function ---
# test_cases = {
#     "Benzene": "c1ccccc1",           # Valid
#     "Ethanol": "CCO",                # Valid
#     "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Valid
#     "Invalid Carbon": "C=C=C=C=",    # Invalid (Valency issue)
#     "Broken Ring": "C1CCC",          # Invalid (Unclosed ring)
#     "Gibberish": "Hello123!"         # Invalid
# }

# print(f"{'Name':<15} | {'SMILES':<15} | {'Valid?'}")
# print("-" * 45)

# for name, smiles in test_cases.items():
#     valid = is_valid_smiles(smiles)
#     print(f"{name:<15} | {smiles:<15} | {valid}")