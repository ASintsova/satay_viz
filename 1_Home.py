import streamlit as st
# import plotly.graph_objects as go
# import pandas as pd
# import plotly.express as px
# import numpy as np
#
# from multipage import MultiPage
# from pages import Analysis, Graphs, Count, Home
#
# app = MultiPage()
# pages = {'Home': ('Home', Home.app),
#          'Count': ('Count', Count.app),
#          'Analysis': ('Analysis', Analysis.app),
#          'Graphs': ('Graphs', Graphs.app)}
#
# for page_name, page in pages.items():
#     app.add_page(page[0], page[1])

st.set_page_config(page_title="SATAY", layout="wide")

#app.run()



def app():
    st.markdown("# SATAY Data Exploration")
    st.markdown("""
                ## Available Pages:

                ### 📊 Intervals
                Compare transposon insertion counts between different conditions using interval-based analysis. 
                Upload differential analysis results to visualize fold changes, significance levels, and generate chromosome-wide scatter plots.

                ### 🔬 Single Replicates  
                Analyze individual replicate datasets with upset plots and ranking visualizations.
                Compare gene classifications across different conditions and explore overlaps between experimental categories.

                ### 🧬 Gene View
                Examine transposon insertions at the individual gene level with detailed genomic context.
                Upload insertion files to visualize insertion patterns around specific genes of interest with customizable padding regions.

                """)
app()
