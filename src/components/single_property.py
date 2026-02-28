import streamlit as st
import pandas as pd

from utils.SMILE_handler import is_valid_smiles, draw_smile_2D, get_smile_features
from utils.main import Models

from .actions import save_single_prediction


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
    """
    )

    smile_input = col2.text_input("Enter SMILES string:", placeholder="e.g., CCO, C1=CC=CC=C1")

    if smile_input:
        if is_valid_smiles(smile_input):
            col2.success("✓ Valid SMILES string")

            img, canon_smiles = draw_smile_2D(smile_input)
            if img:
                col2.image(img, caption=f"Canonical SMILES: {canon_smiles}", width=300)

            try:
                with st.spinner("Predicting properties..."):
                    models = Models()
                    predictions = models.predict_properties(smile_input)
                    features = get_smile_features(smile_input)

                if isinstance(predictions, dict):
                    col2.subheader("Predicted Properties:")

                    prop_cols = col2.columns(len(predictions))
                    for idx, (prop_name, value) in enumerate(predictions.items()):
                        with prop_cols[idx]:
                            st.metric(label=prop_name, value=f"{value:.3f}")

                    col2.subheader("Detailed Results:")
                    df = pd.DataFrame([predictions])
                    col2.dataframe(df, use_container_width=True)

                    st.session_state.current_predictions = predictions
                    st.session_state.current_smiles = canon_smiles
                    st.session_state.current_features = features

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
