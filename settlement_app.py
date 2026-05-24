import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from math import log10

st.set_page_config(page_title="Settlement Prediction Tool", layout="wide")

st.title("Consolidation Settlement Calculator")
st.write("Estimate primary consolidation settlement for normally consolidated clay.")

H = st.number_input("Clay thickness, H (ft)", min_value=0.1, value=30.0)
e0 = st.number_input("Initial void ratio, e₀", min_value=0.0, value=0.65)
sigma0 = st.number_input("Initial effective stress, σ′₀ (psf)", min_value=1.0, value=2910.0)
delta_sigma = st.number_input("Stress increase, Δσ (psf)", min_value=0.0, value=200.0)
Cc = st.number_input("Compression index, Cc", min_value=0.001, value=0.29)

settlement_ft = H * Cc / (1 + e0) * log10((sigma0 + delta_sigma) / sigma0)
settlement_in = settlement_ft * 12

st.subheader("Settlement Result")
st.metric("Settlement", f"{settlement_in:.2f} inches")
st.metric("Settlement", f"{settlement_ft:.4f} ft")

stress_values = np.linspace(50, 1000, 50)
settlement_values = [
    H * Cc / (1 + e0) * log10((sigma0 + ds) / sigma0) * 12
    for ds in stress_values
]

df = pd.DataFrame({
    "Stress Increase Δσ (psf)": stress_values,
    "Settlement (inches)": settlement_values
})

fig = px.line(
    df,
    x="Stress Increase Δσ (psf)",
    y="Settlement (inches)",
    title="Settlement Sensitivity to Stress Increase"
)

st.plotly_chart(fig, use_container_width=True)
st.dataframe(df)
