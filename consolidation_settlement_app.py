import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from math import log10, pi


# ============================================================
# Page Configuration
# ============================================================

st.set_page_config(
    page_title="1D Consolidation Settlement Calculator",
    page_icon="🏗️",
    layout="wide"
)

st.title("1D Consolidation Settlement Calculator with Time History")
st.markdown(
    """
This app calculates **primary consolidation settlement** and estimates the **time required**
to achieve that settlement using classical **Terzaghi 1D consolidation theory**.
"""
)


# ============================================================
# Helper Functions
# ============================================================

def calculate_settlement_nc(H, e0, Cc, sigma_v0, delta_sigma):
    """
    Normally consolidated clay settlement.
    H in ft, stresses in psf.
    Returns settlement in ft.
    """
    sigma_final = sigma_v0 + delta_sigma

    if sigma_v0 <= 0 or sigma_final <= 0:
        return np.nan

    settlement = (H / (1 + e0)) * Cc * log10(sigma_final / sigma_v0)
    return settlement


def calculate_settlement_oc(H, e0, Cc, Cr, sigma_v0, delta_sigma, sigma_p):
    """
    Overconsolidated clay settlement.
    H in ft, stresses in psf.
    Returns settlement in ft.
    """
    sigma_final = sigma_v0 + delta_sigma

    if sigma_v0 <= 0 or sigma_final <= 0 or sigma_p <= 0:
        return np.nan

    if sigma_final <= sigma_p:
        settlement = (H / (1 + e0)) * Cr * log10(sigma_final / sigma_v0)
    else:
        settlement_recompression = (H / (1 + e0)) * Cr * log10(sigma_p / sigma_v0)
        settlement_virgin = (H / (1 + e0)) * Cc * log10(sigma_final / sigma_p)
        settlement = settlement_recompression + settlement_virgin

    return settlement


def average_degree_consolidation_series(Tv, terms=100):
    """
    Average degree of consolidation U using the exact Fourier series solution.

    U = 1 - sum[8 / (pi^2 * (2m+1)^2) * exp(-(2m+1)^2 * pi^2 * Tv / 4)]

    Tv = dimensionless time factor
    Returns U as decimal between 0 and 1.
    """
    if Tv <= 0:
        return 0.0

    summation = 0.0

    for m in range(terms):
        n = 2 * m + 1
        term = (8 / (pi ** 2 * n ** 2)) * np.exp(-(n ** 2 * pi ** 2 * Tv) / 4)
        summation += term

    U = 1 - summation
    return max(0.0, min(1.0, U))


def average_degree_consolidation_approx(Tv):
    """
    Approximate U from Tv.
    This is useful for quick calculations but the app primarily uses the series solution.
    """
    if Tv <= 0:
        return 0.0

    U_values = np.linspace(0.001, 0.999, 1000)
    Tv_values = []

    for U in U_values:
        if U <= 0.60:
            Tv_calc = (pi / 4) * U ** 2
        else:
            Tv_calc = -0.933 * np.log10(1 - U) - 0.085

        Tv_values.append(Tv_calc)

    Tv_values = np.array(Tv_values)

    idx = np.argmin(np.abs(Tv_values - Tv))
    return U_values[idx]


def time_factor_from_U_series(U_target, terms=100):
    """
    Finds Tv corresponding to a target average degree of consolidation U
    using bisection and the series solution.

    U_target should be decimal between 0 and 1.
    """
    if U_target <= 0:
        return 0.0

    if U_target >= 0.9999:
        U_target = 0.9999

    Tv_low = 0.0
    Tv_high = 10.0

    for _ in range(100):
        Tv_mid = 0.5 * (Tv_low + Tv_high)
        U_mid = average_degree_consolidation_series(Tv_mid, terms=terms)

        if U_mid < U_target:
            Tv_low = Tv_mid
        else:
            Tv_high = Tv_mid

    return 0.5 * (Tv_low + Tv_high)


def time_factor_from_U_approx(U):
    """
    Classical approximate relation for time factor Tv from U.
    U should be decimal between 0 and 1.
    """
    if U <= 0:
        return 0.0

    if U < 0.60:
        Tv = (pi / 4) * U ** 2
    else:
        Tv = -0.933 * np.log10(1 - U) - 0.085

    return Tv


def convert_cv_to_ft2_per_day(Cv_input, unit):
    """
    Converts coefficient of consolidation to ft^2/day.
    """
    if unit == "ft²/day":
        return Cv_input

    elif unit == "ft²/year":
        return Cv_input / 365.25

    elif unit == "in²/min":
        return Cv_input * (1 / 144) * 60 * 24

    elif unit == "cm²/sec":
        cm2_to_ft2 = 0.00107639104
        sec_to_day = 86400
        return Cv_input * cm2_to_ft2 * sec_to_day

    elif unit == "m²/year":
        m2_to_ft2 = 10.7639104
        return Cv_input * m2_to_ft2 / 365.25

    else:
        return Cv_input


def convert_time_from_days(days):
    """
    Converts days into a readable set of units.
    """
    years = days / 365.25
    months = days / 30.4375

    return days, months, years


def build_time_history(S_final_ft, Cv_ft2_day, Hdr_ft, max_years, n_points=300):
    """
    Builds time history table for consolidation settlement.
    """
    max_days = max_years * 365.25

    # Log-spaced time avoids crowding early-time values.
    time_days = np.logspace(-3, np.log10(max_days), n_points)

    Tv_values = Cv_ft2_day * time_days / (Hdr_ft ** 2)

    U_values = np.array([
        average_degree_consolidation_series(Tv) for Tv in Tv_values
    ])

    settlement_ft = U_values * S_final_ft
    settlement_in = settlement_ft * 12

    df = pd.DataFrame({
        "Time (days)": time_days,
        "Time (months)": time_days / 30.4375,
        "Time (years)": time_days / 365.25,
        "Time Factor, Tv": Tv_values,
        "Average Degree of Consolidation, U (%)": U_values * 100,
        "Settlement (ft)": settlement_ft,
        "Settlement (in)": settlement_in
    })

    return df


def make_settlement_time_plot(df):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["Time (years)"],
            y=df["Settlement (in)"],
            mode="lines",
            name="Settlement"
        )
    )

    fig.update_layout(
        title="Settlement versus Time",
        xaxis_title="Time (years)",
        yaxis_title="Settlement (in)",
        template="plotly_white",
        height=500
    )

    return fig


def make_settlement_log_time_plot(df):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["Time (days)"],
            y=df["Settlement (in)"],
            mode="lines",
            name="Settlement"
        )
    )

    fig.update_layout(
        title="Settlement versus Time - Log Time Scale",
        xaxis_title="Time (days)",
        yaxis_title="Settlement (in)",
        xaxis_type="log",
        template="plotly_white",
        height=500
    )

    return fig


def make_consolidation_time_plot(df):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["Time (years)"],
            y=df["Average Degree of Consolidation, U (%)"],
            mode="lines",
            name="Average Degree of Consolidation"
        )
    )

    fig.update_layout(
        title="Average Degree of Consolidation versus Time",
        xaxis_title="Time (years)",
        yaxis_title="Average Degree of Consolidation, U (%)",
        template="plotly_white",
        height=500
    )

    return fig


# ============================================================
# Sidebar Inputs
# ============================================================

st.sidebar.header("Input Parameters")

st.sidebar.subheader("Clay Layer Geometry")

H = st.sidebar.number_input(
    "Compressible clay layer thickness, H (ft)",
    min_value=0.1,
    value=30.0,
    step=0.5
)

drainage_condition = st.sidebar.selectbox(
    "Drainage condition",
    [
        "Double drainage - drainage at top and bottom",
        "Single drainage - drainage at one boundary only",
        "Custom drainage path"
    ]
)

if drainage_condition == "Double drainage - drainage at top and bottom":
    Hdr = H / 2
elif drainage_condition == "Single drainage - drainage at one boundary only":
    Hdr = H
else:
    Hdr = st.sidebar.number_input(
        "Drainage path length, Hdr (ft)",
        min_value=0.01,
        value=H / 2,
        step=0.5
    )

st.sidebar.markdown(f"**Drainage path, Hdr = {Hdr:.2f} ft**")


st.sidebar.subheader("Stress Conditions")

stress_input_method = st.sidebar.radio(
    "How do you want to define initial effective stress?",
    [
        "Enter initial vertical effective stress directly",
        "Calculate from effective unit weight and depth"
    ]
)

if stress_input_method == "Enter initial vertical effective stress directly":
    sigma_v0 = st.sidebar.number_input(
        "Initial vertical effective stress, σ'v0 (psf)",
        min_value=1.0,
        value=3000.0,
        step=100.0
    )
else:
    gamma_eff = st.sidebar.number_input(
        "Effective unit weight, γ' (pcf)",
        min_value=1.0,
        value=97.0,
        step=1.0
    )

    z_mid = st.sidebar.number_input(
        "Depth to middle of clay layer, z (ft)",
        min_value=0.1,
        value=30.0,
        step=0.5
    )

    sigma_v0 = gamma_eff * z_mid

delta_sigma = st.sidebar.number_input(
    "Increase in vertical stress, Δσv (psf)",
    min_value=0.0,
    value=200.0,
    step=50.0
)

sigma_final = sigma_v0 + delta_sigma


st.sidebar.subheader("Compressibility Parameters")

soil_condition = st.sidebar.selectbox(
    "Clay consolidation condition",
    [
        "Normally consolidated clay",
        "Overconsolidated clay"
    ]
)

e0 = st.sidebar.number_input(
    "Initial void ratio, e₀",
    min_value=0.01,
    value=0.65,
    step=0.01
)

Cc = st.sidebar.number_input(
    "Compression index, Cc",
    min_value=0.001,
    value=0.29,
    step=0.01
)

if soil_condition == "Overconsolidated clay":
    Cr = st.sidebar.number_input(
        "Recompression index, Cr",
        min_value=0.0001,
        value=0.05,
        step=0.005
    )

    sigma_p = st.sidebar.number_input(
        "Preconsolidation pressure, σ'p (psf)",
        min_value=1.0,
        value=5000.0,
        step=100.0
    )
else:
    Cr = None
    sigma_p = None


st.sidebar.subheader("Time Rate Parameters")

Cv_input = st.sidebar.number_input(
    "Coefficient of consolidation, Cv",
    min_value=0.000001,
    value=1.0,
    step=0.1,
    format="%.6f"
)

Cv_unit = st.sidebar.selectbox(
    "Cv unit",
    [
        "ft²/day",
        "ft²/year",
        "in²/min",
        "cm²/sec",
        "m²/year"
    ]
)

Cv_ft2_day = convert_cv_to_ft2_per_day(Cv_input, Cv_unit)

max_years = st.sidebar.number_input(
    "Maximum time for time-history plot (years)",
    min_value=0.01,
    value=30.0,
    step=1.0
)

target_U_percent = st.sidebar.slider(
    "Target average degree of consolidation, U (%)",
    min_value=1,
    max_value=99,
    value=90,
    step=1
)

series_terms = st.sidebar.slider(
    "Number of Fourier series terms",
    min_value=10,
    max_value=300,
    value=100,
    step=10
)


# ============================================================
# Settlement Calculation
# ============================================================

if soil_condition == "Normally consolidated clay":
    S_final_ft = calculate_settlement_nc(
        H=H,
        e0=e0,
        Cc=Cc,
        sigma_v0=sigma_v0,
        delta_sigma=delta_sigma
    )
else:
    S_final_ft = calculate_settlement_oc(
        H=H,
        e0=e0,
        Cc=Cc,
        Cr=Cr,
        sigma_v0=sigma_v0,
        delta_sigma=delta_sigma,
        sigma_p=sigma_p
    )

S_final_in = S_final_ft * 12


# ============================================================
# Time History Calculation
# ============================================================

df_time = build_time_history(
    S_final_ft=S_final_ft,
    Cv_ft2_day=Cv_ft2_day,
    Hdr_ft=Hdr,
    max_years=max_years,
    n_points=350
)

target_U_decimal = target_U_percent / 100

Tv_target_series = time_factor_from_U_series(
    target_U_decimal,
    terms=series_terms
)

Tv_target_approx = time_factor_from_U_approx(target_U_decimal)

time_target_days_series = Tv_target_series * (Hdr ** 2) / Cv_ft2_day
time_target_days_approx = Tv_target_approx * (Hdr ** 2) / Cv_ft2_day

target_settlement_ft = target_U_decimal * S_final_ft
target_settlement_in = target_settlement_ft * 12

days_series, months_series, years_series = convert_time_from_days(time_target_days_series)
days_approx, months_approx, years_approx = convert_time_from_days(time_target_days_approx)


# ============================================================
# Main Results
# ============================================================

st.header("Calculation Results")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "Ultimate Primary Settlement",
        f"{S_final_in:.3f} in"
    )

with col2:
    st.metric(
        "Initial Effective Stress, σ'v0",
        f"{sigma_v0:,.0f} psf"
    )

with col3:
    st.metric(
        "Final Effective Stress, σ'vf",
        f"{sigma_final:,.0f} psf"
    )

with col4:
    st.metric(
        "Cv Used",
        f"{Cv_ft2_day:.6f} ft²/day"
    )


st.subheader("Time to Reach Target Consolidation")

col5, col6, col7, col8 = st.columns(4)

with col5:
    st.metric(
        "Target U",
        f"{target_U_percent:.0f}%"
    )

with col6:
    st.metric(
        "Settlement at Target U",
        f"{target_settlement_in:.3f} in"
    )

with col7:
    st.metric(
        "Time Required",
        f"{years_series:.2f} years"
    )

with col8:
    st.metric(
        "Time Required",
        f"{days_series:,.0f} days"
    )


st.info(
    f"""
Using the Fourier series solution, the estimated time to reach **{target_U_percent:.0f}%**
average consolidation is approximately **{days_series:,.0f} days**, or **{years_series:.2f} years**.
At that time, the estimated settlement is **{target_settlement_in:.3f} inches**.
"""
)


# ============================================================
# Optional Target Settlement Calculation
# ============================================================

st.header("Time to Reach a User-Selected Settlement")

target_settlement_user_in = st.number_input(
    "Enter target settlement amount (in)",
    min_value=0.0,
    max_value=max(float(S_final_in), 0.001),
    value=float(0.9 * S_final_in) if S_final_in > 0 else 0.0,
    step=0.1
)

if S_final_in > 0:
    U_from_settlement = target_settlement_user_in / S_final_in
    U_from_settlement = min(max(U_from_settlement, 0.0), 0.9999)

    Tv_from_settlement = time_factor_from_U_series(
        U_from_settlement,
        terms=series_terms
    )

    time_from_settlement_days = Tv_from_settlement * (Hdr ** 2) / Cv_ft2_day
    d_set, m_set, y_set = convert_time_from_days(time_from_settlement_days)

    col9, col10, col11 = st.columns(3)

    with col9:
        st.metric(
            "Equivalent U",
            f"{U_from_settlement * 100:.1f}%"
        )

    with col10:
        st.metric(
            "Time",
            f"{d_set:,.0f} days"
        )

    with col11:
        st.metric(
            "Time",
            f"{y_set:.2f} years"
        )

else:
    st.warning("Ultimate settlement is zero or invalid. Check the input parameters.")


# ============================================================
# Plots
# ============================================================

st.header("Time-History Plots")

tab1, tab2, tab3 = st.tabs(
    [
        "Settlement vs Time",
        "Settlement vs Log Time",
        "Degree of Consolidation vs Time"
    ]
)

with tab1:
    fig_settlement = make_settlement_time_plot(df_time)

    fig_settlement.add_trace(
        go.Scatter(
            x=[years_series],
            y=[target_settlement_in],
            mode="markers+text",
            text=[f"{target_U_percent:.0f}% U"],
            textposition="top center",
            name="Target Point"
        )
    )

    st.plotly_chart(fig_settlement, use_container_width=True)

with tab2:
    fig_log = make_settlement_log_time_plot(df_time)

    fig_log.add_trace(
        go.Scatter(
            x=[days_series],
            y=[target_settlement_in],
            mode="markers+text",
            text=[f"{target_U_percent:.0f}% U"],
            textposition="top center",
            name="Target Point"
        )
    )

    st.plotly_chart(fig_log, use_container_width=True)

with tab3:
    fig_U = make_consolidation_time_plot(df_time)

    fig_U.add_trace(
        go.Scatter(
            x=[years_series],
            y=[target_U_percent],
            mode="markers+text",
            text=[f"{target_U_percent:.0f}% U"],
            textposition="top center",
            name="Target Point"
        )
    )

    st.plotly_chart(fig_U, use_container_width=True)


# ============================================================
# Time History Table
# ============================================================

st.header("Time-History Table")

display_df = df_time.copy()

display_df["Time (days)"] = display_df["Time (days)"].round(3)
display_df["Time (months)"] = display_df["Time (months)"].round(3)
display_df["Time (years)"] = display_df["Time (years)"].round(3)
display_df["Time Factor, Tv"] = display_df["Time Factor, Tv"].round(5)
display_df["Average Degree of Consolidation, U (%)"] = display_df[
    "Average Degree of Consolidation, U (%)"
].round(2)
display_df["Settlement (ft)"] = display_df["Settlement (ft)"].round(5)
display_df["Settlement (in)"] = display_df["Settlement (in)"].round(4)

st.dataframe(display_df, use_container_width=True)

csv = display_df.to_csv(index=False).encode("utf-8")

st.download_button(
    label="Download Time-History Table as CSV",
    data=csv,
    file_name="consolidation_time_history.csv",
    mime="text/csv"
)


# ============================================================
# Engineering Interpretation
# ============================================================

st.header("Engineering Interpretation")

if S_final_in < 0:
    st.error("Computed settlement is negative. Please check stress and soil inputs.")
elif S_final_in == 0:
    st.warning("Computed settlement is zero. This may occur if stress increase is zero.")
else:
    st.markdown(
        f"""
### Summary

- Clay layer thickness, **H** = `{H:.2f} ft`
- Drainage path, **Hdr** = `{Hdr:.2f} ft`
- Initial effective vertical stress, **σ'v0** = `{sigma_v0:,.0f} psf`
- Stress increase, **Δσv** = `{delta_sigma:,.0f} psf`
- Final effective vertical stress, **σ'vf** = `{sigma_final:,.0f} psf`
- Ultimate primary consolidation settlement = **{S_final_in:.3f} in**
- Estimated time to reach **{target_U_percent:.0f}% consolidation** = **{years_series:.2f} years**
- Settlement at **{target_U_percent:.0f}% consolidation** = **{target_settlement_in:.3f} in**

### Practical Note

The calculated primary consolidation settlement is the long-term theoretical settlement
associated with dissipation of excess pore water pressure. In practice, the observed
settlement rate can vary depending on soil layering, drainage boundaries, secondary
compression, construction staging, preload/surcharge conditions, and the actual field value
of **Cv**.
"""
    )


# ============================================================
# Formula Reference
# ============================================================

st.header("Formula Reference")

with st.expander("Primary Consolidation Settlement Formulas", expanded=True):

    st.markdown(
        r"""
### Normally Consolidated Clay

\[
S_c = \frac{H}{1 + e_0} C_c \log_{10}
\left(
\frac{\sigma'_{v0} + \Delta \sigma_v}{\sigma'_{v0}}
\right)
\]

Where:

- \( S_c \) = primary consolidation settlement  
- \( H \) = thickness of compressible clay layer  
- \( e_0 \) = initial void ratio  
- \( C_c \) = compression index  
- \( \sigma'_{v0} \) = initial vertical effective stress  
- \( \Delta \sigma_v \) = increase in vertical stress  

---

### Overconsolidated Clay

If final stress is less than or equal to preconsolidation pressure:

\[
S_c = \frac{H}{1 + e_0} C_r \log_{10}
\left(
\frac{\sigma'_{vf}}{\sigma'_{v0}}
\right)
\]

If final stress exceeds preconsolidation pressure:

\[
S_c =
\frac{H}{1 + e_0} C_r \log_{10}
\left(
\frac{\sigma'_p}{\sigma'_{v0}}
\right)
+
\frac{H}{1 + e_0} C_c \log_{10}
\left(
\frac{\sigma'_{vf}}{\sigma'_p}
\right)
\]

Where:

\[
\sigma'_{vf} = \sigma'_{v0} + \Delta \sigma_v
\]

- \( C_r \) = recompression index  
- \( \sigma'_p \) = preconsolidation pressure  
"""
    )

with st.expander("Time Rate of Consolidation Formulas", expanded=True):

    st.markdown(
        r"""
### Dimensionless Time Factor

\[
T_v = \frac{C_v t}{H_{dr}^2}
\]

Solving for time:

\[
t = \frac{T_v H_{dr}^2}{C_v}
\]

Where:

- \( T_v \) = dimensionless time factor  
- \( C_v \) = coefficient of consolidation  
- \( t \) = elapsed time  
- \( H_{dr} \) = maximum drainage path length  

---

### Drainage Path

For single drainage:

\[
H_{dr} = H
\]

For double drainage:

\[
H_{dr} = \frac{H}{2}
\]

---

### Average Degree of Consolidation

The app uses the Fourier series solution:

\[
U = 1 -
\sum_{m=0}^{\infty}
\frac{8}{\pi^2(2m+1)^2}
\exp
\left[
-\frac{(2m+1)^2 \pi^2 T_v}{4}
\right]
\]

Where:

- \( U \) = average degree of consolidation as a decimal  
- \( T_v \) = time factor  

---

### Settlement at Time t

\[
S(t) = U(t) S_c
\]

Where:

- \( S(t) \) = settlement at time \( t \)  
- \( U(t) \) = average degree of consolidation at time \( t \)  
- \( S_c \) = ultimate primary consolidation settlement  
"""
    )


# ============================================================
# Footer
# ============================================================

st.markdown("---")
st.caption(
    "This app is intended for preliminary geotechnical engineering calculations. "
    "Final design should consider project-specific subsurface data, laboratory consolidation test results, "
    "field instrumentation, construction staging, and applicable geotechnical standards of practice."
)
