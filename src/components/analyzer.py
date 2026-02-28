import streamlit as st
import pandas as pd

from utils.SMILE_handler import is_valid_smiles
from utils.analysis import Analyser

from .actions import save_analyzer_results


def render_analyzer(col2):
    """Render Analyzer for Multiple SMILES."""
    col2.header("Analyzer - Properties and Plots for Multiple SMILES")

    col2.markdown("""
    Enter multiple SMILES strings to analyze and visualize their properties.
    The analyzer will:
    - Predict properties for all valid SMILES
    - Sort them by priority order
    - Generate visualization plots
    """
    )

    input_method = col2.radio("Input Method:", ["Text Input", "CSV Upload"])
    smiles_list = []

    if input_method == "Text Input":
        smiles_text = col2.text_area(
            "Enter SMILES strings (one per line):",
            placeholder="CCO\nC1=CC=CC=C1\nCCC\n...",
        )
        if smiles_text:
            smiles_list = [s.strip() for s in smiles_text.split("\n") if s.strip()]
    else:
        csv_file = col2.file_uploader("Upload CSV file (must have 'SMILES' column):", type=["csv"])
        if csv_file:
            df = pd.read_csv(csv_file)
            if "SMILES" in df.columns:
                smiles_list = df["SMILES"].tolist()
            else:
                col2.error("CSV must contain a 'SMILES' column")

    if smiles_list:
        col2.info(f"Found {len(smiles_list)} SMILES strings")

        valid_smiles = [s for s in smiles_list if is_valid_smiles(s)]
        col2.write(f"Valid SMILES: {len(valid_smiles)} / {len(smiles_list)}")

        if valid_smiles:
            col2.subheader("Select Priority Order:")
            all_properties = ["Tg", "Tc", "Rg", "FFV", "Density"]
            priority_order = col2.multiselect(
                "Drag to reorder (select in priority order):",
                options=all_properties,
                default=all_properties,
            )

            if priority_order:
                try:
                    with st.spinner("Analyzing SMILES..."):
                        analyser = Analyser(valid_smiles)
                        plots, sorted_results = analyser.get(priority_order)

                    col2.subheader("Analysis Results:")
                    col2.write(f"Analyzed {len(sorted_results)} compounds")

                    col2.subheader("Sorted SMILES with Predicted Properties:")
                    if sorted_results:
                        results_df = pd.DataFrame(sorted_results)
                        cols = results_df.columns.tolist()
                        if "SMILES" in cols:
                            cols.remove("SMILES")
                            results_df = results_df[["SMILES"] + cols]
                        results_df.insert(0, "Rank", range(1, len(results_df) + 1))
                        for col in ["Tg", "Tc", "Rg", "FFV", "Density"]:
                            if col in results_df.columns:
                                results_df[col] = results_df[col].apply(lambda x: f"{x:.4f}")
                        col2.dataframe(results_df, use_container_width=True)

                        st.session_state.analyzer_results = sorted_results

                        col2.divider()
                        save_col1, save_col2 = col2.columns([1, 3])
                        with save_col1:
                            if st.button("💾 Save All Results to Database", key="save_analyzer_results"):
                                save_analyzer_results(col2, sorted_results)

                    if plots:
                        col2.subheader("Property Distribution Plots:")
                        for prop_name, fig in plots.items():
                            col2.pyplot(fig)

                except Exception as e:
                    col2.error(f"Error during analysis: {str(e)}")
        else:
            col2.warning("No valid SMILES found in the input.")
