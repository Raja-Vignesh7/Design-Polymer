import streamlit as st

from utils.main import db_config_info
from utils.database import DatabaseHandler


def render_settings(col2):
    """Render Settings for Database Configuration."""
    col2.header("Settings - Database Configuration")

    col2.markdown("""
    Configure your database connection parameters here.
    These settings will be saved to the .env file.
    """
    )

    tab1, tab2 = col2.tabs(["Configure Database", "Connection Test"])

    with tab1:
        st.subheader("Database Connection Parameters")
        col_a, col_b = st.columns(2)
        with col_a:
            host = st.text_input("Database Host:", value="localhost", placeholder="localhost or IP address")
            user = st.text_input("Username:", value="root", placeholder="MySQL username")
            password = st.text_input("Password:", type="password", placeholder="MySQL password")
        with col_b:
            db_name = st.text_input("Database Name:", placeholder="Database name")
            port = st.number_input("Port:", value=3306, min_value=1, max_value=65535)

        if st.button("Save Configuration"):
            try:
                config = db_config_info()
                config.set_config_info(host, user, password, db_name, port)
                st.session_state.db_config = {
                    "DB_HOST": host,
                    "DB_USER": user,
                    "DB_PASSWORD": password,
                    "DB_NAME": db_name,
                    "DB_PORT": port,
                }
                st.success("✓ Configuration saved successfully!")
            except Exception as e:
                st.error(f"Error saving configuration: {str(e)}")

    with tab2:
        st.subheader("Test Database Connection")

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
                        port=int(current_config.get("DB_PORT", 3306)),
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
