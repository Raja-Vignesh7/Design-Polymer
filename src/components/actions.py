import streamlit as st

from utils.SMILE_handler import get_smile_features
from utils.main import db_config_info
from utils.database import DatabaseHandler


def save_single_prediction(container, smiles, predictions, features):
    """Save a single polymer prediction to the database.

    Args:
        container: Streamlit container for displaying messages
        smiles (str): Canonical SMILES string
        predictions (dict): Dictionary of predicted properties
        features (array-like): Molecular features
    """
    try:
        config_manager = db_config_info()
        config = config_manager.get_config_info()

        handler = DatabaseHandler(
            host=config.get("DB_HOST", "localhost"),
            user=config.get("DB_USER", "root"),
            password=config.get("DB_PASSWORD", ""),
            database=config.get("DB_NAME", "polymer_db"),
            port=int(config.get("DB_PORT", 3306)),
        )

        if not handler.connect():
            container.error("❌ Failed to connect to database. Please check your configuration in Settings.")
            return

        success = handler.save_polymer_data(
            smiles=smiles,
            features=features,
            tg=float(predictions.get("Tg", 0)),
            ffv=float(predictions.get("FFV", 0)),
            density=float(predictions.get("Density", 0)),
            tc=float(predictions.get("Tc", 0)),
            rg=float(predictions.get("Rg", 0)),
        )

        handler.disconnect()

        if success:
            container.success(f"✓ Successfully saved polymer data to database!\n\nSMILES: {smiles}")
        else:
            container.error("❌ Failed to save data to database.")
    except Exception as e:
        container.error(f"❌ Error saving to database: {str(e)}")


def save_analyzer_results(container, results):
    """Save multiple polymer predictions to the database.

    Args:
        container: Streamlit container for displaying messages
        results (list): List of dictionaries containing SMILES and predicted properties
    """
    try:
        config_manager = db_config_info()
        config = config_manager.get_config_info()
        handler = DatabaseHandler(
            host=config.get("DB_HOST", "localhost"),
            user=config.get("DB_USER", "root"),
            password=config.get("DB_PASSWORD", ""),
            database=config.get("DB_NAME", "polymer_db"),
            port=int(config.get("DB_PORT", 3306)),
        )

        if not handler.connect():
            container.error("❌ Failed to connect to database. Please check your configuration in Settings.")
            return

        saved_count = 0
        failed_count = 0

        with st.spinner("Saving results to database..."):
            for result in results:
                try:
                    smiles = result.get("SMILES", "")
                    features = get_smile_features(smiles, as_dict=True)
                    success = handler.save_polymer_data(
                        smiles=smiles,
                        features=features,
                        tg=float(result.get("Tg", 0)),
                        ffv=float(result.get("FFV", 0)),
                        density=float(result.get("Density", 0)),
                        tc=float(result.get("Tc", 0)),
                        rg=float(result.get("Rg", 0)),
                    )
                    if success:
                        saved_count += 1
                    else:
                        failed_count += 1
                except Exception as row_error:
                    container.warning(f"⚠ Failed to save record for SMILES: {smiles}. Error: {str(row_error)}")
                    failed_count += 1

        handler.disconnect()

        if saved_count > 0:
            container.success(f"✓ Successfully saved {saved_count} polymer records to database!")
        if failed_count > 0:
            container.warning(f"⚠ Failed to save {failed_count} records.")
    except Exception as e:
        container.error(f"❌ Error saving to database: {str(e)}")
