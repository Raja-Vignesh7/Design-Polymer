import streamlit as st
from utils.SMILE_handler import is_valid_smiles, draw_smile_2D

def main():
    st.set_page_config(page_title="Design Polymer", page_icon=":alembic:", layout="wide")

    col1, col2, col3 = st.columns([1, 6, 1])
    
    with col2:
        st.title("Design Polymer :alembic:")
        smile_iput = st.text_input("Enter a SMILES string for a polymer fragment:", "")
        if smile_iput:
            if is_valid_smiles(smile_iput):
                img, canon_smiles = draw_smile_2D(smile_iput)
                st.image(img, caption=f"Canonical SMILES: {canon_smiles}")
            else:
                st.error("Invalid SMILES string. Please check your input.")
        st.markdown(
            """
            Welcome to the Design Polymer application! This tool allows you to design and analyze polymers with ease.
            
            **Features:**
            - Interactive polymer design interface
            - Real-time property calculations
            - Visualization of polymer structures
            
            **Getting Started:**
            1. Use the sidebar to navigate through different sections.
            2. Input your desired polymer parameters.
            3. Analyze and visualize the results.
            
            Enjoy designing your polymers!
            """
        )
    with col1:
        # st.metric(label="Version", value="1.0.0")
        st.header("Menu")
        
        
    with col3:
        st.header("Info")
        st.markdown(
            """
            **Design Polymer App**
            
            - **Version:** 1.0.0"""
        )
        
if __name__ == "__main__":
    main()