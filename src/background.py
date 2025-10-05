"""
Background component for serving Next.js app as Streamlit background.
Designed for Azure deployment.
"""
import streamlit as st


def render_nextjs_background():
    """
    Embeds the exact Next.js WebGL shader background via iframe.
    Requires Next.js running on port 3000.
    """

    # Inject the iframe and styles directly into the page
    st.markdown("""
    <iframe
        id="nextjs-background-iframe"
        src="http://localhost:3000"
        style="
            position: fixed;
            top: 0;
            left: 0;
            width: 100vw;
            height: 100vh;
            border: none;
            z-index: -9999;
            pointer-events: none;
        "
    ></iframe>

    <style>
        /* Make Streamlit app transparent to show iframe behind */
        .stApp {
            background: transparent !important;
        }

        /* Ensure main container is above background */
        .main {
            position: relative;
            z-index: 1;
        }

        /* Semi-transparent main content area */
        .main .block-container {
            background: rgba(0, 0, 0, 0.5) !important;
            backdrop-filter: blur(15px) !important;
            -webkit-backdrop-filter: blur(15px) !important;
            border-radius: 16px !important;
            padding: 2rem !important;
            margin-top: 2rem !important;
        }

        /* Semi-transparent sidebar */
        [data-testid="stSidebar"] {
            background: rgba(0, 0, 0, 0.75) !important;
            backdrop-filter: blur(12px) !important;
            -webkit-backdrop-filter: blur(12px) !important;
            position: relative;
            z-index: 2;
        }

        [data-testid="stSidebar"] > div:first-child {
            background: transparent !important;
        }

        /* Text visibility */
        h1, h2, h3, h4, h5, h6 {
            color: white !important;
            text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.5);
        }

        .stMarkdown, p, span, label {
            color: white !important;
        }

        /* Metrics styling */
        .stMetric {
            background: rgba(255, 255, 255, 0.15) !important;
            backdrop-filter: blur(10px) !important;
            padding: 1rem !important;
            border-radius: 8px !important;
            border: 1px solid rgba(255, 255, 255, 0.2) !important;
        }

        .stMetric label {
            color: rgba(255, 255, 255, 0.9) !important;
        }

        .stMetric [data-testid="stMetricValue"] {
            color: white !important;
        }

        /* Tabs */
        .stTabs [data-baseweb="tab-list"] {
            background: rgba(0, 0, 0, 0.3) !important;
            border-radius: 8px !important;
        }

        .stTabs [data-baseweb="tab"] {
            color: rgba(255, 255, 255, 0.8) !important;
        }

        .stTabs [aria-selected="true"] {
            color: white !important;
            background: rgba(59, 130, 246, 0.4) !important;
        }

        /* Inputs */
        input, textarea, select {
            background: rgba(255, 255, 255, 0.1) !important;
            color: white !important;
            border: 1px solid rgba(255, 255, 255, 0.3) !important;
        }

        /* Buttons */
        .stButton > button {
            background: rgba(59, 130, 246, 0.6) !important;
            color: white !important;
            border: 1px solid rgba(255, 255, 255, 0.3) !important;
        }

        .stButton > button:hover {
            background: rgba(59, 130, 246, 0.8) !important;
        }

        /* Ensure header is above background */
        header {
            position: relative;
            z-index: 3;
        }
    </style>
    """, unsafe_allow_html=True)
