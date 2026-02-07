import streamlit as st
import pandas as pd
import sys
import os

# Add utils to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'utils'))

from utils.SMILE_handler import is_valid_smiles, draw_smile_2D, get_smile_features
from utils.main import Models, db_config_info
from utils.analysis import Analyser
from utils.database import DatabaseHandler

def save_single_prediction(container, smiles, predictions, features):
    """
    Save a single polymer prediction to the database.
    
    Args:
        container: Streamlit container for displaying messages
        smiles (str): Canonical SMILES string
        predictions (dict): Dictionary of predicted properties
        features (array-like): Molecular features
    """
    try:
        # Load database configuration
        config_manager = db_config_info()
        config = config_manager.get_config_info()
        
        # Initialize database handler
        handler = DatabaseHandler(
            host=config.get("DB_HOST", "localhost"),
            user=config.get("DB_USER", "root"),
            password=config.get("DB_PASSWORD", ""),
            database=config.get("DB_NAME", "polymer_db"),
            port=int(config.get("DB_PORT", 3306))
        )
        
        # Connect to database
        if not handler.connect():
            container.error("❌ Failed to connect to database. Please check your configuration in Settings.")
            return
        
        # Save the polymer data
        success = handler.save_polymer_data(
            smiles=smiles,
            features=features,
            tg=float(predictions.get("Tg", 0)),
            ffv=float(predictions.get("FFV", 0)),
            density=float(predictions.get("Density", 0)),
            tc=float(predictions.get("Tc", 0)),
            rg=float(predictions.get("Rg", 0))
        )
        
        # Disconnect from database
        handler.disconnect()
        
        if success:
            container.success(f"✓ Successfully saved polymer data to database!\n\nSMILES: {smiles}")
        else:
            container.error("❌ Failed to save data to database.")
    
    except Exception as e:
        container.error(f"❌ Error saving to database: {str(e)}")

def save_analyzer_results(container, results):
    """
    Save multiple polymer predictions to the database.
    
    Args:
        container: Streamlit container for displaying messages
        results (list): List of dictionaries containing SMILES and predicted properties
    """
    try:
        # Load database configuration
        config_manager = db_config_info()
        config = config_manager.get_config_info()
        
        # Initialize database handler
        handler = DatabaseHandler(
            host=config.get("DB_HOST", "localhost"),
            user=config.get("DB_USER", "root"),
            password=config.get("DB_PASSWORD", ""),
            database=config.get("DB_NAME", "polymer_db"),
            port=int(config.get("DB_PORT", 3306))
        )
        
        # Connect to database
        if not handler.connect():
            container.error("❌ Failed to connect to database. Please check your configuration in Settings.")
            return
        
        # Save each result to the database
        saved_count = 0
        failed_count = 0
        
        with st.spinner("Saving results to database..."):
            for result in results:
                try:
                    smiles = result.get("SMILES", "")
                    
                    # Extract features for this SMILES if not already present
                    if "features" in result:
                        features = result.get("features", {})
                    else:
                        features = get_smile_features(smiles)
                    
                    success = handler.save_polymer_data(
                        smiles=smiles,
                        features=features,
                        tg=float(result.get("Tg", 0)),
                        ffv=float(result.get("FFV", 0)),
                        density=float(result.get("Density", 0)),
                        tc=float(result.get("Tc", 0)),
                        rg=float(result.get("Rg", 0))
                    )
                    
                    if success:
                        saved_count += 1
                    else:
                        failed_count += 1
                except Exception as row_error:
                    failed_count += 1
        
        # Disconnect from database
        handler.disconnect()
        
        # Display results
        if saved_count > 0:
            container.success(f"✓ Successfully saved {saved_count} polymer records to database!")
        
        if failed_count > 0:
            container.warning(f"⚠ Failed to save {failed_count} records.")

    except Exception as e:
        container.error(f"❌ Error saving to database: {str(e)}")

def initialize_session_state():
    """Initialize Streamlit session state variables."""
    if 'current_menu' not in st.session_state:
        st.session_state.current_menu = "Get Properties"
    if 'db_config' not in st.session_state:
        st.session_state.db_config = None

def render_menu(col1):
    """Render the menu options in column 1."""
    st.session_state.current_menu = col1.radio(
        "Select Menu Option:",
        options=[
            "Get Properties",
            "Analyzer",
            "Settings"
        ],
        index=0
    )

def render_get_properties(col2):
    """Render Get Properties for One Formula interface."""
    col2.header("Get Properties for One Formula")
    
    col2.markdown("""
    Enter a SMILES string to predict polymer properties including:
    - **Tg**: Glass Transition Temperature
    - **Tc**: Thermal Conductivity
    - **Rg**: Radius of Gyration
    - **FFV**: Free Volume Fraction
    - **Density**: Material Density
    """)
    
    smile_input = col2.text_input("Enter SMILES string:", placeholder="e.g., CCO, C1=CC=CC=C1")
    
    if smile_input:
        if is_valid_smiles(smile_input):
            col2.success("✓ Valid SMILES string")
            
            # Display molecular structure
            img, canon_smiles = draw_smile_2D(smile_input)
            if img:
                col2.image(img, caption=f"Canonical SMILES: {canon_smiles}", width=300)
            
            # Get predictions
            try:
                with st.spinner("Predicting properties..."):
                    models = Models()
                    predictions = models.predict_properties(smile_input)
                    features = get_smile_features(smile_input)
                
                if isinstance(predictions, dict):
                    col2.subheader("Predicted Properties:")
                    
                    # Display as metrics
                    prop_cols = col2.columns(len(predictions))
                    for idx, (prop_name, value) in enumerate(predictions.items()):
                        with prop_cols[idx]:
                            st.metric(label=prop_name, value=f"{value:.3f}")
                    
                    # Display as table
                    col2.subheader("Detailed Results:")
                    df = pd.DataFrame([predictions])
                    col2.dataframe(df, use_container_width=True)
                    
                    # Store predictions in session state
                    st.session_state.current_predictions = predictions
                    st.session_state.current_smiles = canon_smiles
                    st.session_state.current_features = features
                    
                    # Add save to database button
                    col2.divider()
                    save_col1, save_col2 = col2.columns([1, 3])
                    
                    with save_col1:
                        if st.button("💾 Save to Database", key="save_single_prediction"):
                            save_single_prediction(col2, canon_smiles, predictions, features)
                    
            except Exception as e:
                col2.error(f"Error predicting properties: {str(e)}")
        else:
            col2.error("❌ Invalid SMILES string. Please check your input.")
    else:
        col2.info("Enter a SMILES string above to see predictions.")

def render_analyzer(col2):
    """Render Analyzer for Multiple SMILES."""
    col2.header("Analyzer - Properties and Plots for Multiple SMILES")
    
    col2.markdown("""
    Enter multiple SMILES strings to analyze and visualize their properties.
    The analyzer will:
    - Predict properties for all valid SMILES
    - Sort them by priority order
    - Generate visualization plots
    """)
    
    # Input methods
    input_method = col2.radio("Input Method:", ["Text Input", "CSV Upload"])
    
    smiles_list = []
    
    if input_method == "Text Input":
        smiles_text = col2.text_area(
            "Enter SMILES strings (one per line):",
            placeholder="CCO\nC1=CC=CC=C1\nCCC\n..."
        )
        if smiles_text:
            smiles_list = [s.strip() for s in smiles_text.split('\n') if s.strip()]
    else:
        csv_file = col2.file_uploader("Upload CSV file (must have 'SMILES' column):", type=['csv'])
        if csv_file:
            df = pd.read_csv(csv_file)
            if 'SMILES' in df.columns:
                smiles_list = df['SMILES'].tolist()
            else:
                col2.error("CSV must contain a 'SMILES' column")
    
    if smiles_list:
        col2.info(f"Found {len(smiles_list)} SMILES strings")
        
        # Validate SMILES
        valid_smiles = [s for s in smiles_list if is_valid_smiles(s)]
        col2.write(f"Valid SMILES: {len(valid_smiles)} / {len(smiles_list)}")
        
        if valid_smiles:
            # Priority order selection
            col2.subheader("Select Priority Order:")
            all_properties = ["Tg", "Tc", "Rg", "FFV", "Density"]
            priority_order = col2.multiselect(
                "Drag to reorder (select in priority order):",
                options=all_properties,
                default=all_properties
            )
            
            if priority_order:
                try:
                    with st.spinner("Analyzing SMILES..."):
                        analyser = Analyser(valid_smiles)
                        plots, sorted_results = analyser.get(priority_order)
                    
                    col2.subheader("Analysis Results:")
                    col2.write(f"Analyzed {len(sorted_results)} compounds")
                    
                    # Display sorted SMILES with their properties
                    col2.subheader("Sorted SMILES with Predicted Properties:")
                    if sorted_results:
                        results_df = pd.DataFrame(sorted_results)
                        # Reorder columns to put SMILES first
                        cols = results_df.columns.tolist()
                        if 'SMILES' in cols:
                            cols.remove('SMILES')
                            results_df = results_df[['SMILES'] + cols]
                        
                        # Add rank column
                        results_df.insert(0, 'Rank', range(1, len(results_df) + 1))
                        
                        # Format numeric columns to 4 decimal places
                        for col in ['Tg', 'Tc', 'Rg', 'FFV', 'Density']:
                            if col in results_df.columns:
                                results_df[col] = results_df[col].apply(lambda x: f"{x:.4f}")
                        
                        col2.dataframe(results_df, use_container_width=True)
                        
                        # Store results in session state
                        st.session_state.analyzer_results = sorted_results
                        
                        # Add save to database button
                        col2.divider()
                        save_col1, save_col2 = col2.columns([1, 3])
                        
                        with save_col1:
                            if st.button("💾 Save All Results to Database", key="save_analyzer_results"):
                                save_analyzer_results(col2, sorted_results)
                    
                    # Display plots
                    if plots:
                        col2.subheader("Property Distribution Plots:")
                        for prop_name, fig in plots.items():
                            col2.pyplot(fig)
                    
                except Exception as e:
                    col2.error(f"Error during analysis: {str(e)}")
        else:
            col2.warning("No valid SMILES found in the input.")

def render_settings(col2):
    """Render Settings for Database Configuration."""
    col2.header("Settings - Database Configuration")
    
    col2.markdown("""
    Configure your database connection parameters here.
    These settings will be saved to the .env file.
    """)
    
    # Create tabs for database settings
    tab1, tab2 = col2.tabs(["Configure Database", "Connection Test"])
    
    with tab1:
        st.subheader("Database Connection Parameters")
        
        col_a, col_b = st.columns(2)
        
        with col_a:
            host = st.text_input("Database Host:", value="localhost", placeholder="localhost or IP address")
            user = st.text_input("Username:", value="root", placeholder="MySQL username")
            password = st.text_input("Password:", type="password", placeholder="MySQL password")
        
        with col_b:
            db_name = st.text_input("Database Name:", value="polymer_db", placeholder="Database name")
            port = st.number_input("Port:", value=3306, min_value=1, max_value=65535)
        
        # Save button
        if st.button("Save Configuration"):
            try:
                config = db_config_info()
                config.set_config_info(host, user, password, db_name, port)
                st.session_state.db_config = {
                    "DB_HOST": host,
                    "DB_USER": user,
                    "DB_PASSWORD": password,
                    "DB_NAME": db_name,
                    "DB_PORT": port
                }
                st.success("✓ Configuration saved successfully!")
            except Exception as e:
                st.error(f"Error saving configuration: {str(e)}")
    
    with tab2:
        st.subheader("Test Database Connection")
        
        # Load current configuration
        try:
            config_manager = db_config_info()
            current_config = config_manager.get_config_info()
            
            st.info("Testing with current configuration...")
            
            if st.button("Test Connection"):
                try:
                    handler = DatabaseHandler(
                        host=current_config.get("DB_HOST", "localhost"),
                        user=current_config.get("DB_USER", "root"),
                        password=current_config.get("DB_PASSWORD", ""),
                        database=current_config.get("DB_NAME", "polymer_db"),
                        port=int(current_config.get("DB_PORT", 3306))
                    )
                    
                    if handler.connect():
                        st.success("✓ Connection successful!")
                        handler.disconnect()
                    else:
                        st.error("✗ Connection failed. Check your credentials.")
                except Exception as e:
                    st.error(f"Error during connection test: {str(e)}")
        except Exception as e:
            st.warning(f"No configuration found. Please configure the database first. ({str(e)})")

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
    """)
    
    col3.divider()
    col3.subheader("Menu Guide")
    col3.markdown("""
    **Get Properties**
    Submit a single SMILES to predict its properties
    
    **Analyzer**
    Analyze multiple SMILES and visualize results
    
    **Settings**
    Configure database connection
    """)

def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Design Polymer",
        page_icon=":alembic:",
        layout="wide",
        initial_sidebar_state="collapsed"
    )
    
    # Initialize session state
    initialize_session_state()
    
    # Main title
    st.title("Design Polymer :alembic:")
    
    # Create 3-column layout
    col1, col2, col3 = st.columns([1.2, 5, 1.5])
    
    # Column 1: Menu
    with col1:
        st.subheader("Menu")
        render_menu(col1)
    
    # Column 2: Content
    with col2:
        if st.session_state.current_menu == "Get Properties":
            render_get_properties(col2)
        elif st.session_state.current_menu == "Analyzer":
            render_analyzer(col2)
        elif st.session_state.current_menu == "Settings":
            render_settings(col2)
    
    # Column 3: About
    with col3:
        render_about(col3)

if __name__ == "__main__":
    main()