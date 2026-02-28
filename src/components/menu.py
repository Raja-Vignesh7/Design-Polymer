import streamlit as st


def initialize_session_state():
    """Ensure the Streamlit session state has our expected keys."""
    if "current_menu" not in st.session_state:
        st.session_state.current_menu = "Get Properties"
    if "db_config" not in st.session_state:
        st.session_state.db_config = None


def render_menu(col1):
    """Render the navigation menu in the left column."""
    st.session_state.current_menu = col1.selectbox(
        "Select Menu Option:",
        options=[
            "Get Properties",
            "Analyzer",
            "Settings",
        ],
        index=0,
    )
