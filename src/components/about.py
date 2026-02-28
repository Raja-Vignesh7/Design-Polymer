import streamlit as st


def render_about(col3):
    """Render About/Info section in column 3."""
    col3.header("About")
    col3.markdown("""
    ### Design Polymer
    
    **Version:** 1.0.0
    
    #### Features
    - 🧬 Single Formula Analysis
    - 📊 Multiple SMILES Analysis
    - 📈 Property Predictions
    - 🗄️ Database Integration
    - ⚙️ Configuration Management
    
    #### Capabilities
    The application predicts five key polymer properties:
    - **Tg**: Glass Transition Temperature
    - **Tc**: Thermal Conductivity  
    - **Rg**: Radius of Gyration
    - **FFV**: Free Volume Fraction
    - **Density**: Material Density
    
    #### Supported Input
    - SMILES strings
    - CSV files
    - Direct text input
    
    #### Technologies
    - Streamlit (UI)
    - RDKit (Chemistry)
    - Mordred (Descriptors)
    - scikit-learn (ML)
    - MySQL (Database)
    """
    )

    col3.divider()
    col3.subheader("Menu Guide")
    col3.markdown("""
    **Get Properties**
    Submit a single SMILES to predict its properties
    
    **Analyzer**
    Analyze multiple SMILES and visualize results
    
    **Settings**
    Configure database connection
    """
    )
