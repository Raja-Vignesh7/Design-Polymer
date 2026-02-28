import streamlit as st
import pandas as pd
import sys
import os

# ensure the utils package is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "utils"))

# ui components
from components.menu import initialize_session_state, render_menu
from components.single_property import render_get_properties
from components.analyzer import render_analyzer
from components.settings import render_settings
from components.about import render_about


def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Design Polymer",
        page_icon=":alembic:",
        layout="wide",
        initial_sidebar_state="collapsed",
    )

    # make sure session state keys exist
    initialize_session_state()

    # header
    st.title("Design Polymer :alembic:")

    # three-column layout
    col1, col2, col3 = st.columns([1.2, 5, 1.5])

    with col1:
        st.subheader("Menu")
        render_menu(col1)

    with col2:
        if st.session_state.current_menu == "Get Properties":
            render_get_properties(col2)
        elif st.session_state.current_menu == "Analyzer":
            render_analyzer(col2)
        elif st.session_state.current_menu == "Settings":
            render_settings(col2)

    with col3:
        render_about(col3)


if __name__ == "__main__":
    main()
