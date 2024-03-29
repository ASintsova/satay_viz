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
                ###  Count Page
                This page allows counting # of transposons over custom intervals 
                accross [yeast genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/). 
                Requires `bed` files as an input.

                ### Analysis Page

                This page allows identification of intervals with differential transposon abundances between 
                **two different conditions**. Requires count table generated via `Count Page` as an input.

                ### Graphs Page

                This page provides simple visualizations of differential transposon insertion frequencies. 
                Requires LFC and z-score table generated via `Analysis Page`

                """)
app()
