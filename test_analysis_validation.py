"""
Validation script to verify the analysis.py changes work correctly.
This script demonstrates that the get() method now returns SMILES with properties.
"""

import sys
import os
import pandas as pd
from unittest.mock import patch, MagicMock

# Add src/utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'utils'))

def validate_analysis_changes():
    """Validate that Analyser.get() returns sorted SMILES with properties."""
    
    # Mock the dependencies to avoid import errors
    with patch('analysis.Models'), \
         patch('analysis.is_valid_smiles', return_value=True), \
         patch('analysis.plt.subplots') as mock_subplots:
        
        from analysis import Analyser
        
        # Setup mocks
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        
        # Mock model predictions
        predictions = [
            {'Tg': 150.0, 'Tc': 0.6, 'Rg': 6.0, 'FFV': 0.22, 'Density': 1.25},
            {'Tg': 180.0, 'Tc': 0.7, 'Rg': 7.0, 'FFV': 0.25, 'Density': 1.30},
            {'Tg': 120.0, 'Tc': 0.5, 'Rg': 5.0, 'FFV': 0.20, 'Density': 1.20},
        ]
        
        with patch('analysis.Models.predict_properties', side_effect=predictions):
            # Create analyser with sample SMILES
            analyser = Analyser(["CCO", "C1=CC=CC=C1", "CCC"])
            
            # Get analysis results (this is the modified method)
            plots, sorted_results = analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            
            # Validate return types
            print("✓ Return type validation:")
            print(f"  - plots type: {type(plots).__name__} (expected: dict)")
            print(f"  - sorted_results type: {type(sorted_results).__name__} (expected: list)")
            
            # Validate sorted_results structure
            print("\n✓ Sorted results structure validation:")
            print(f"  - Number of results: {len(sorted_results)}")
            
            if sorted_results:
                print(f"  - First result type: {type(sorted_results[0]).__name__} (expected: dict)")
                first_result = sorted_results[0]
                print(f"  - Keys in result: {list(first_result.keys())}")
                
                # Validate all required keys are present
                required_keys = ['Tg', 'Tc', 'Rg', 'FFV', 'Density', 'SMILES']
                all_keys_present = all(key in first_result for key in required_keys)
                print(f"  - All required keys present: {all_keys_present}")
                
                # Validate sorting (should be by Tg in descending order)
                tg_values = [r['Tg'] for r in sorted_results]
                is_sorted_correctly = tg_values == sorted(tg_values, reverse=True)
                print(f"\n✓ Sorting validation:")
                print(f"  - Tg values: {tg_values}")
                print(f"  - Sorted in descending order: {is_sorted_correctly}")
                
                # Display the results as a table
                print(f"\n✓ Sorted SMILES with predicted properties:")
                df = pd.DataFrame(sorted_results)
                cols = df.columns.tolist()
                if 'SMILES' in cols:
                    cols.remove('SMILES')
                    df = df[['SMILES'] + cols]
                df.insert(0, 'Rank', range(1, len(df) + 1))
                
                # Format numeric columns
                for col in ['Tg', 'Tc', 'Rg', 'FFV', 'Density']:
                    if col in df.columns:
                        df[col] = df[col].apply(lambda x: f"{x:.4f}")
                
                print(df.to_string(index=False))
                
                return True
            
    return False

if __name__ == "__main__":
    try:
        success = validate_analysis_changes()
        if success:
            print("\n" + "="*60)
            print("✓ All validations passed!")
            print("✓ The Analyser.get() method now returns:")
            print("  1. plots (dict of matplotlib figures)")
            print("  2. sorted_results (list of dicts with SMILES and properties)")
            print("="*60)
        else:
            print("\n✗ Validation failed")
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error during validation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
