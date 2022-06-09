import streamlit as st

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