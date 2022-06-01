import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import numpy as np

from multipage import MultiPage
from pages import Analysis, Graphs, Count

app = MultiPage()
pages = {'Count': ('Count', Count.app),
         'Analysis': ('Analysis', Analysis.app),
         'Graphs': ('Graphs', Graphs.app)}

for page_name, page in pages.items():
    app.add_page(page[0], page[1])

st.set_page_config(page_title="SATAY", layout="wide")

app.run()