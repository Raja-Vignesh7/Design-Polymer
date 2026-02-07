import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from SMILE_handler import is_valid_smiles, draw_smile_2D
from main import Models

class Analyser:
    def __init__(self, data):
        """
        Initialize the Analyser with SMILES data.
        
        Args:
            data (str or list): Single SMILES string or list of SMILES strings
        """
        # Normalize input to always work with list internally
        if isinstance(data, str):
            self.smiles_list = [data]
        else:
            self.smiles_list = data if isinstance(data, list) else [data]
        
        # Initialize models for predictions
        self.models = Models()
        
        # Store predictions and SMILES together
        self.predictions = []
        self.valid_smiles = []
        
        # Get predictions for all SMILES
        for smi in self.smiles_list:
            if is_valid_smiles(smi):
                pred = self.models.predict_properties(smi)
                self.predictions.append(pred)
                self.valid_smiles.append(smi)

    def get(self, priority_order):
        """
        Analyze SMILES data and return sorted results with visualizations.
        
        Creates bar plots for each target property and returns SMILES sorted by priority order.
        
        Args:
            priority_order (list): List of property names in priority order.
                                  First element is highest priority, last is lowest.
                                  Valid properties: ["Tg", "Tc", "Rg", "FFV", "Density"]
                                  
        Example:
            >>> analyser = Analyser(["CCO", "C1=CC=CC=C1", "CCC"])
            >>> plots, sorted_smiles = analyser.get(["Tg", "Tc", "Rg", "FFV", "Density"])
            >>> # plots is a dict of matplotlib figures
            >>> # sorted_smiles is a list of SMILES sorted by Tg (descending)
            
        Returns:
            tuple: (plots_dict, sorted_smiles)
                - plots_dict (dict): Dictionary with property names as keys and 
                                    matplotlib Figure objects as values
                - sorted_smiles (list): List of SMILES sorted by priority order
                
        Notes:
            - Sorts by the first property in priority_order (primary sort key)
            - Higher values are considered better (descending sort)
            - Only valid SMILES are included in results
        """
        if not self.valid_smiles or not self.predictions:
            print("Error: No valid SMILES to analyze.")
            return {}, []
        
        # Create a DataFrame for easier manipulation
        df = pd.DataFrame(self.predictions)
        df['SMILES'] = self.valid_smiles
        
        # Sort by the primary priority (first element in priority_order)
        for property in priority_order[::-1]:
            df = df.sort_values(by=property, ascending=False).reset_index(drop=True)
        
        # Get sorted SMILES list
        sorted_smiles = df['SMILES'].tolist()
        
        # Create bar plots for each target property
        plots = {}
        valid_properties = ["Tg", "Tc", "Rg", "FFV", "Density"]
        
        # Set style for better-looking plots
        sns.set_style("whitegrid")
        
        for prop in valid_properties:
            if prop in df.columns:
                # Create figure
                fig, ax = plt.subplots(figsize=(12, 6))
                
                # Create bar plot
                bars = ax.bar(range(len(df)), df[prop], color='steelblue', alpha=0.8, edgecolor='navy')
                
                # Add value labels on top of bars
                for i, (bar, val) in enumerate(zip(bars, df[prop])):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                           f'{val:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
                
                # Customize plot
                ax.set_xlabel('SMILES Compounds', fontsize=12, fontweight='bold')
                ax.set_ylabel(f'{prop} Value', fontsize=12, fontweight='bold')
                ax.set_title(f'{prop} Distribution Across SMILES', fontsize=14, fontweight='bold')
                ax.set_xticks(range(len(df)))
                ax.set_xticklabels([f'SMILES {i+1}' for i in range(len(df))], rotation=45, ha='right')
                
                # Highlight primary property with different color
                
                for bar in bars:
                    bar.set_color('coral')
                    bar.set_alpha(0.9)
                ax.set_title(f'{prop} Distribution (PRIMARY PROPERTY)', fontsize=14, fontweight='bold', color='red')
                
                plt.tight_layout()
                plots[prop] = fig
        
        print(f"\n✓ Analysis Complete!")
        print(f"✓ Generated {len(plots)} plot(s) for properties: {', '.join(plots.keys())}")
        print(f"✓ SMILES sorted by primary priority: {priority_order}")
        print(f"✓ Total valid SMILES analyzed: {len(sorted_smiles)}")
        
        return plots, sorted_smiles
        