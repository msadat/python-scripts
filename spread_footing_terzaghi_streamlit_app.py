import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# ============================================================
# Streamlit Page Setup
# ============================================================
st.set_page_config(
    page_title="Spread Footing Design - Terzaghi Bearing Capacity",
    page_icon="🏗️",
    layout="wide",
)


# ============================================================
# Engineering Helper Functions
# ============================================================
@dataclass
class BearingCapacityFactors:
    nc: float
    nq: float
    ngamma: float


@dataclass
class BearingCapacityResult:
    qult_gross_ksf: float
    qult_net_ksf: float
    qall_gross_ksf: float
    qall_net_ksf: float
    surcharge_q_ksf: float
    factors: BearingCapacityFactors
    shape_c: float
    shape_q: float
    shape_gamma: float


@dataclass
class PressureResult:
    q_avg_ksf: float
    q_min_ksf: float
    q_max_ksf: float
    e_x_ft: float
    e_y_ft: float
    kern_ok_x: bool
    kern_ok_y: bool
    no_tension: bool
    uplift_area_present: bool
    grid_x: np.ndarray
    grid_y: np.ndarray
    grid_q: np.ndarray


def safe_float(value: float, default: float = 0.0) -> float:
    try:
        if value is None or math.isnan(float(value)):
            return default
        return float(value)
    except Exception:
        return default


def bearing_capacity_factors(phi_deg: float, ngamma_method: str) -> BearingCapacityFactors:
    """
    Bearing capacity factors using classical closed-form expressions.

    For phi = 0 undrained condition, Terzaghi's commonly used values are:
        Nc = 5.7, Nq = 1.0, Ngamma = 0.0

    For phi > 0:
        Nq = exp(pi tan phi) tan^2(45 + phi/2)
        Nc = (Nq - 1) / tan phi

    Ngamma is not unique across classical methods, so the app lets the user select:
        - Terzaghi approximation: Ngamma = 1.5 (Nq - 1) tan phi
        - Meyerhof:             Ngamma = (Nq - 1) tan(1.4 phi)
        - Vesic:                Ngamma = 2 (Nq + 1) tan phi
    """
    phi_deg = max(0.0, min(phi_deg, 50.0))

    if abs(phi_deg) < 1e-9:
        return BearingCapacityFactors(nc=5.7, nq=1.0, ngamma=0.0)

    phi = math.radians(phi_deg)
    nq = math.exp(math.pi * math.tan(phi)) * math.tan(math.radians(45.0) + phi / 2.0) ** 2
    nc = (nq - 1.0) / math.tan(phi)

    if ngamma_method == "Meyerhof: (Nq - 1) tan(1.4φ)":
        ngamma = (nq - 1.0) * math.tan(1.4 * phi)
    elif ngamma_method == "Vesic: 2(Nq + 1) tanφ":
        ngamma = 2.0 * (nq + 1.0) * math.tan(phi)
    else:
        ngamma = 1.5 * (nq - 1.0) * math.tan(phi)

    return BearingCapacityFactors(nc=nc, nq=nq, ngamma=ngamma)


def effective_unit_weight_for_bearing(
    gamma_moist_pcf: float,
    gamma_sat_pcf: float,
    water_table_depth_ft: float,
    foundation_depth_ft: float,
    footing_width_ft: float,
) -> Tuple[float, str]:
    """
    Simplified groundwater correction for the gamma term.

    Depth is measured from ground surface downward.
    - If water table is below Df + B, use moist unit weight.
    - If water table is at or above Df, use submerged unit weight for bearing gamma term.
    - If water table lies between Df and Df + B, interpolate between moist and submerged.

    This correction affects the 0.5 gamma B Ngamma term. Surcharge q is handled separately.
    """
    gamma_w = 62.4
    gamma_sub_pcf = max(gamma_sat_pcf - gamma_w, 0.0)

    if water_table_depth_ft <= foundation_depth_ft:
        return gamma_sub_pcf, "Water table at/above foundation base; submerged unit weight used for bearing γ term."

    if water_table_depth_ft >= foundation_depth_ft + footing_width_ft:
        return gamma_moist_pcf, "Water table below Df + B; moist unit weight used for bearing γ term."

    distance_below_base = water_table_depth_ft - foundation_depth_ft
    ratio_moist = max(0.0, min(distance_below_base / footing_width_ft, 1.0))
    gamma_eff = ratio_moist * gamma_moist_pcf + (1.0 - ratio_moist) * gamma_sub_pcf
    return gamma_eff, "Water table between Df and Df + B; interpolated effective γ used for bearing γ term."


def surcharge_at_base(
    gamma_moist_pcf: float,
    gamma_sat_pcf: float,
    water_table_depth_ft: float,
    foundation_depth_ft: float,
) -> float:
    """
    Computes approximate effective overburden surcharge q at foundation base in ksf.
    """
    gamma_w = 62.4
    gamma_sub_pcf = max(gamma_sat_pcf - gamma_w, 0.0)

    if foundation_depth_ft <= 0:
        return 0.0

    if water_table_depth_ft >= foundation_depth_ft:
        q_psf = gamma_moist_pcf * foundation_depth_ft
    elif water_table_depth_ft <= 0:
        q_psf = gamma_sub_pcf * foundation_depth_ft
    else:
        q_psf = gamma_moist_pcf * water_table_depth_ft + gamma_sub_pcf * (foundation_depth_ft - water_table_depth_ft)

    return q_psf / 1000.0


def shape_factors_terzaghi_rectangular(B_ft: float, L_ft: float) -> Tuple[float, float, float]:
    """
    Terzaghi-style shape modifiers for rectangular/square footings.

    For rectangular footing with B <= L:
        sc = 1 + 0.3 B/L
        sq = 1.0
        sγ = 1 - 0.2 B/L

    This recovers approximately:
        square: sc = 1.3, sγ = 0.8
    Terzaghi square expression is commonly written as:
        qult = 1.3 cNc + qNq + 0.4γBNγ
    which is equivalent to 0.5γBNγ * 0.8.
    """
    B = max(min(B_ft, L_ft), 1e-9)
    L = max(max(B_ft, L_ft), 1e-9)
    ratio = B / L
    sc = 1.0 + 0.3 * ratio
    sq = 1.0
    sgamma = 1.0 - 0.2 * ratio
    return sc, sq, sgamma


def compute_bearing_capacity(
    c_psf: float,
    phi_deg: float,
    gamma_moist_pcf: float,
    gamma_sat_pcf: float,
    water_table_depth_ft: float,
    Df_ft: float,
    B_ft: float,
    L_ft: float,
    fs_target: float,
    ngamma_method: str,
) -> BearingCapacityResult:
    B_eff = max(min(B_ft, L_ft), 1e-9)
    factors = bearing_capacity_factors(phi_deg, ngamma_method)
    sc, sq, sgamma = shape_factors_terzaghi_rectangular(B_ft, L_ft)

    q_ksf = surcharge_at_base(gamma_moist_pcf, gamma_sat_pcf, water_table_depth_ft, Df_ft)
    gamma_eff_pcf, _ = effective_unit_weight_for_bearing(
        gamma_moist_pcf, gamma_sat_pcf, water_table_depth_ft, Df_ft, B_eff
    )

    c_ksf = c_psf / 1000.0
    gamma_eff_kcf = gamma_eff_pcf / 1000.0

    qult_gross = (
        c_ksf * factors.nc * sc
        + q_ksf * factors.nq * sq
        + 0.5 * gamma_eff_kcf * B_eff * factors.ngamma * sgamma
    )
    qult_net = max(qult_gross - q_ksf, 0.0)
    qall_gross = qult_gross / fs_target if fs_target > 0 else 0.0
    qall_net = qult_net / fs_target if fs_target > 0 else 0.0

    return BearingCapacityResult(
        qult_gross_ksf=qult_gross,
        qult_net_ksf=qult_net,
        qall_gross_ksf=qall_gross,
        qall_net_ksf=qall_net,
        surcharge_q_ksf=q_ksf,
        factors=factors,
        shape_c=sc,
        shape_q=sq,
        shape_gamma=sgamma,
    )


def compute_pressure_distribution(
    P_kips: float,
    Mx_kipft: float,
    My_kipft: float,
    B_ft: float,
    L_ft: float,
    grid_n: int = 81,
) -> PressureResult:
    """
    Rigid footing linear contact pressure under axial load and biaxial moment.

    Coordinate system:
        x = footing width direction, from -B/2 to +B/2
        y = footing length direction, from -L/2 to +L/2

    Moments:
        Mx causes pressure variation along y.
        My causes pressure variation along x.

    q(x,y) = P/A + My*x/Iy + Mx*y/Ix

    where:
        Ix = B L^3 / 12
        Iy = L B^3 / 12
    """
    B = max(B_ft, 1e-9)
    L = max(L_ft, 1e-9)
    A = B * L
    Ix = B * L**3 / 12.0
    Iy = L * B**3 / 12.0

    x = np.linspace(-B / 2.0, B / 2.0, grid_n)
    y = np.linspace(-L / 2.0, L / 2.0, grid_n)
    X, Y = np.meshgrid(x, y)

    q_avg = P_kips / A if A > 0 else 0.0
    q_grid = q_avg + (My_kipft * X / Iy if Iy > 0 else 0.0) + (Mx_kipft * Y / Ix if Ix > 0 else 0.0)

    e_x = My_kipft / P_kips if abs(P_kips) > 1e-9 else 0.0
    e_y = Mx_kipft / P_kips if abs(P_kips) > 1e-9 else 0.0

    kern_ok_x = abs(e_x) <= B / 6.0
    kern_ok_y = abs(e_y) <= L / 6.0
    q_min = float(np.min(q_grid))
    q_max = float(np.max(q_grid))
    uplift_area_present = q_min < 0.0

    return PressureResult(
        q_avg_ksf=q_avg,
        q_min_ksf=q_min,
        q_max_ksf=q_max,
        e_x_ft=e_x,
        e_y_ft=e_y,
        kern_ok_x=kern_ok_x,
        kern_ok_y=kern_ok_y,
        no_tension=(kern_ok_x and kern_ok_y and not uplift_area_present),
        uplift_area_present=uplift_area_present,
        grid_x=X,
        grid_y=Y,
        grid_q=q_grid,
    )


def build_pressure_surface_figure(pr: PressureResult, B_ft: float, L_ft: float) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Surface(
            x=pr.grid_x,
            y=pr.grid_y,
            z=pr.grid_q,
            colorscale="Viridis",
            colorbar=dict(title="q (ksf)"),
        )
    )
    fig.update_layout(
        title="3D Bearing Pressure Distribution Under Footing",
        scene=dict(
            xaxis_title="Width direction, x (ft)",
            yaxis_title="Length direction, y (ft)",
            zaxis_title="Bearing pressure, q (ksf)",
        ),
        margin=dict(l=0, r=0, b=0, t=45),
        height=620,
    )
    return fig


def build_pressure_contour_figure(pr: PressureResult) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=pr.grid_x[0, :],
            y=pr.grid_y[:, 0],
            z=pr.grid_q,
            colorscale="Viridis",
            colorbar=dict(title="q (ksf)"),
            contours=dict(showlabels=True, labelfont=dict(size=11)),
        )
    )
    fig.update_layout(
        title="Plan View Pressure Contours",
        xaxis_title="Width direction, x (ft)",
        yaxis_title="Length direction, y (ft)",
        height=560,
    )
    return fig


def classify_status(fs_actual: float, fs_target: float, qmax: float, qall: float, uplift: bool) -> Tuple[str, str]:
    if uplift:
        return (
            "⚠️ Not adequate for classical full-contact pressure assumption",
            "Negative bearing pressure is predicted at part of the footing. Soil cannot take tension; increase footing size, reduce moment, or perform partial-contact/eccentric footing analysis.",
        )
    if qmax <= qall and fs_actual >= fs_target:
        return (
            "✅ Adequate for bearing pressure based on selected criteria",
            "Maximum service bearing pressure is less than or equal to allowable bearing pressure.",
        )
    return (
        "❌ Not adequate for bearing capacity",
        "Maximum service bearing pressure exceeds allowable bearing pressure. Increase B/L, improve soil, increase embedment, or reduce loads.",
    )


# ============================================================
# App Header
# ============================================================
st.title("Spread Footing Design Using Terzaghi Bearing Capacity Theory")
st.caption(
    "Preliminary geotechnical design tool for rectangular spread footings under axial load and biaxial moments. "
    "Use for screening and educational calculations; final design should be reviewed by a qualified geotechnical/structural engineer."
)


# ============================================================
# Sidebar Inputs
# ============================================================
st.sidebar.header("1. Loads")
P_kips = st.sidebar.number_input("Service axial load, P (kips)", min_value=0.01, value=250.0, step=10.0)
Mx_kipft = st.sidebar.number_input("Service moment about x-axis, Mx (kip-ft)", value=100.0, step=10.0)
My_kipft = st.sidebar.number_input("Service moment about y-axis, My (kip-ft)", value=75.0, step=10.0)

st.sidebar.header("2. Footing Geometry")
B_ft = st.sidebar.number_input("Footing width, B (ft)", min_value=0.5, value=8.0, step=0.5)
L_ft = st.sidebar.number_input("Footing length, L (ft)", min_value=0.5, value=10.0, step=0.5)
Df_ft = st.sidebar.number_input("Foundation embedment depth, Df (ft)", min_value=0.0, value=3.0, step=0.5)

st.sidebar.header("3. Soil Parameters")
soil_type = st.sidebar.selectbox(
    "Soil behavior / design condition",
    [
        "Drained sand / granular soil",
        "Drained clay or c-φ soil",
        "Undrained clay, φ = 0",
    ],
)

if soil_type == "Undrained clay, φ = 0":
    phi_deg = 0.0
    c_psf = st.sidebar.number_input("Undrained shear strength, su = c (psf)", min_value=0.0, value=1500.0, step=100.0)
elif soil_type == "Drained sand / granular soil":
    phi_deg = st.sidebar.number_input("Effective friction angle, φ' (degrees)", min_value=0.0, max_value=50.0, value=32.0, step=1.0)
    c_psf = st.sidebar.number_input("Effective cohesion, c' (psf)", min_value=0.0, value=0.0, step=50.0)
else:
    phi_deg = st.sidebar.number_input("Effective friction angle, φ' (degrees)", min_value=0.0, max_value=50.0, value=28.0, step=1.0)
    c_psf = st.sidebar.number_input("Effective cohesion, c' (psf)", min_value=0.0, value=500.0, step=50.0)

gamma_moist_pcf = st.sidebar.number_input("Moist unit weight, γmoist (pcf)", min_value=50.0, value=120.0, step=1.0)
gamma_sat_pcf = st.sidebar.number_input("Saturated unit weight, γsat (pcf)", min_value=50.0, value=125.0, step=1.0)
water_table_depth_ft = st.sidebar.number_input(
    "Depth to groundwater below ground surface (ft)", min_value=0.0, value=999.0, step=1.0
)

st.sidebar.header("4. Design Criteria")
fs_target = st.sidebar.number_input("Target factor of safety", min_value=1.0, value=3.0, step=0.1)
capacity_basis = st.sidebar.radio(
    "Bearing comparison basis",
    ["Gross allowable vs gross service qmax", "Net allowable vs net service qmax"],
)
ngamma_method = st.sidebar.selectbox(
    "Nγ method",
    [
        "Terzaghi approximation: 1.5(Nq - 1)tanφ",
        "Meyerhof: (Nq - 1) tan(1.4φ)",
        "Vesic: 2(Nq + 1) tanφ",
    ],
)

st.sidebar.header("5. Plot Settings")
grid_n = st.sidebar.slider("Pressure plot resolution", min_value=25, max_value=151, value=81, step=2)


# ============================================================
# Calculations
# ============================================================
B_eff_for_bearing = min(B_ft, L_ft)
A_ft2 = B_ft * L_ft

bc = compute_bearing_capacity(
    c_psf=c_psf,
    phi_deg=phi_deg,
    gamma_moist_pcf=gamma_moist_pcf,
    gamma_sat_pcf=gamma_sat_pcf,
    water_table_depth_ft=water_table_depth_ft,
    Df_ft=Df_ft,
    B_ft=B_ft,
    L_ft=L_ft,
    fs_target=fs_target,
    ngamma_method=ngamma_method,
)

pr = compute_pressure_distribution(
    P_kips=P_kips,
    Mx_kipft=Mx_kipft,
    My_kipft=My_kipft,
    B_ft=B_ft,
    L_ft=L_ft,
    grid_n=grid_n,
)

gamma_eff_pcf, gw_note = effective_unit_weight_for_bearing(
    gamma_moist_pcf, gamma_sat_pcf, water_table_depth_ft, Df_ft, B_eff_for_bearing
)

service_net_qmax = max(pr.q_max_ksf - bc.surcharge_q_ksf, 0.0)
service_net_qavg = max(pr.q_avg_ksf - bc.surcharge_q_ksf, 0.0)

if capacity_basis == "Gross allowable vs gross service qmax":
    demand_ksf = pr.q_max_ksf
    allowable_ksf = bc.qall_gross_ksf
    ultimate_ksf = bc.qult_gross_ksf
    basis_label = "gross"
else:
    demand_ksf = service_net_qmax
    allowable_ksf = bc.qall_net_ksf
    ultimate_ksf = bc.qult_net_ksf
    basis_label = "net"

fs_actual = ultimate_ksf / demand_ksf if demand_ksf > 1e-9 else float("inf")
status_title, status_body = classify_status(fs_actual, fs_target, demand_ksf, allowable_ksf, pr.uplift_area_present)


# ============================================================
# Main Results
# ============================================================
st.subheader("Design Summary")

summary_col1, summary_col2, summary_col3, summary_col4 = st.columns(4)
summary_col1.metric("Footing Area", f"{A_ft2:,.1f} ft²")
summary_col2.metric("Average Service Pressure", f"{pr.q_avg_ksf:,.2f} ksf")
summary_col3.metric("Maximum Service Pressure", f"{pr.q_max_ksf:,.2f} ksf")
summary_col4.metric("Minimum Service Pressure", f"{pr.q_min_ksf:,.2f} ksf")

status_box = st.container(border=True)
with status_box:
    st.markdown(f"### {status_title}")
    st.write(status_body)
    st.write(
        f"**Selected comparison:** {basis_label.capitalize()} demand = `{demand_ksf:,.2f} ksf`; "
        f"{basis_label.capitalize()} allowable = `{allowable_ksf:,.2f} ksf`; "
        f"Calculated FS = `{fs_actual:,.2f}` vs target FS = `{fs_target:,.2f}`."
    )

if not pr.kern_ok_x or not pr.kern_ok_y:
    st.warning(
        "The resultant is outside the footing kern in at least one direction. "
        "The full-contact linear pressure equation may predict tension. Consider increasing footing size, reducing moments, or using partial-contact analysis."
    )

if pr.uplift_area_present:
    st.error(
        "Negative contact pressure is present in the calculated pressure grid. Soil cannot resist tension. "
        "Do not use the full-contact qmin/qmax result as a final design check."
    )


# ============================================================
# Detailed Output Tables
# ============================================================
st.subheader("Bearing Capacity Calculation")

calc_data = {
    "Parameter": [
        "Nc",
        "Nq",
        "Nγ",
        "Shape factor sc",
        "Shape factor sq",
        "Shape factor sγ",
        "Effective surcharge at base, q",
        "Effective γ for bearing term",
        "Ultimate gross bearing capacity, qult,gross",
        "Ultimate net bearing capacity, qult,net",
        "Allowable gross bearing pressure, qall,gross",
        "Allowable net bearing pressure, qall,net",
    ],
    "Value": [
        f"{bc.factors.nc:,.3f}",
        f"{bc.factors.nq:,.3f}",
        f"{bc.factors.ngamma:,.3f}",
        f"{bc.shape_c:,.3f}",
        f"{bc.shape_q:,.3f}",
        f"{bc.shape_gamma:,.3f}",
        f"{bc.surcharge_q_ksf:,.3f} ksf",
        f"{gamma_eff_pcf:,.1f} pcf",
        f"{bc.qult_gross_ksf:,.3f} ksf",
        f"{bc.qult_net_ksf:,.3f} ksf",
        f"{bc.qall_gross_ksf:,.3f} ksf",
        f"{bc.qall_net_ksf:,.3f} ksf",
    ],
}
st.dataframe(pd.DataFrame(calc_data), use_container_width=True, hide_index=True)
st.info(gw_note)

st.subheader("Eccentricity and Contact Pressure Check")

pressure_data = {
    "Parameter": [
        "e_x = My / P",
        "e_y = Mx / P",
        "B / 6 kern limit",
        "L / 6 kern limit",
        "Average gross pressure",
        "Maximum gross pressure",
        "Minimum gross pressure",
        "Maximum net pressure used in net check",
    ],
    "Value": [
        f"{pr.e_x_ft:,.3f} ft",
        f"{pr.e_y_ft:,.3f} ft",
        f"{B_ft / 6.0:,.3f} ft",
        f"{L_ft / 6.0:,.3f} ft",
        f"{pr.q_avg_ksf:,.3f} ksf",
        f"{pr.q_max_ksf:,.3f} ksf",
        f"{pr.q_min_ksf:,.3f} ksf",
        f"{service_net_qmax:,.3f} ksf",
    ],
    "Status": [
        "OK" if pr.kern_ok_x else "Outside kern",
        "OK" if pr.kern_ok_y else "Outside kern",
        "Reference",
        "Reference",
        "Reference",
        "Demand",
        "No tension OK" if pr.q_min_ksf >= 0 else "Tension predicted",
        "Demand" if capacity_basis.startswith("Net") else "Reference",
    ],
}
st.dataframe(pd.DataFrame(pressure_data), use_container_width=True, hide_index=True)


# ============================================================
# Plots
# ============================================================
st.subheader("Pressure Distribution Under Foundation")
plot_col1, plot_col2 = st.columns(2)

with plot_col1:
    st.plotly_chart(build_pressure_contour_figure(pr), use_container_width=True)

with plot_col2:
    st.plotly_chart(build_pressure_surface_figure(pr, B_ft, L_ft), use_container_width=True)


# ============================================================
# Design Interpretation
# ============================================================
st.subheader("Engineering Interpretation")

interpretation = []
interpretation.append(
    f"The resultant eccentricity is **ex = {pr.e_x_ft:.3f} ft** across the footing width and **ey = {pr.e_y_ft:.3f} ft** along the footing length."
)

if pr.kern_ok_x and pr.kern_ok_y:
    interpretation.append(
        "The resultant lies within the rectangular kern limits in both directions, so the full-contact linear pressure assumption is acceptable for preliminary checks."
    )
else:
    interpretation.append(
        "The resultant is outside at least one kern limit. This indicates possible partial contact and invalidates a simple full-contact linear pressure check for final design."
    )

if demand_ksf <= allowable_ksf and not pr.uplift_area_present:
    interpretation.append(
        "The selected footing size satisfies the selected bearing capacity check. Settlement, sliding, overturning, structural shear, punching shear, and flexure still need to be checked separately."
    )
else:
    interpretation.append(
        "The selected footing size does not satisfy the selected bearing capacity/contact pressure check. Increasing B and/or L is usually the first practical iteration."
    )

for item in interpretation:
    st.markdown(f"- {item}")


# ============================================================
# Formula Reference Section
# ============================================================
st.divider()
st.header("Formula Reference")

with st.expander("1. Terzaghi Bearing Capacity Equation", expanded=True):
    st.latex(r"q_{ult,gross}=cN_c s_c + qN_q s_q + \frac{1}{2}\gamma B N_\gamma s_\gamma")
    st.markdown(
        "Where `B` is the smaller footing dimension for the bearing capacity term, "
        "`q = γDf` is the effective overburden surcharge at the foundation base, "
        "and `sc`, `sq`, and `sγ` are shape factors."
    )
    st.latex(r"q_{ult,net}=q_{ult,gross}-q")
    st.latex(r"q_{all}=\frac{q_{ult}}{FS}")

with st.expander("2. Bearing Capacity Factors"):
    st.latex(r"N_q=e^{\pi \tan\phi}\tan^2\left(45^\circ+\frac{\phi}{2}\right)")
    st.latex(r"N_c=\frac{N_q-1}{\tan\phi}")
    st.markdown("For φ = 0 undrained clay, this app uses:")
    st.latex(r"N_c=5.7,\quad N_q=1.0,\quad N_\gamma=0.0")
    st.markdown("Selectable Nγ options:")
    st.latex(r"N_\gamma=1.5(N_q-1)\tan\phi\quad\text{Terzaghi approximation}")
    st.latex(r"N_\gamma=(N_q-1)\tan(1.4\phi)\quad\text{Meyerhof}")
    st.latex(r"N_\gamma=2(N_q+1)\tan\phi\quad\text{Vesic}")

with st.expander("3. Rectangular Footing Shape Factors"):
    st.markdown("For rectangular footing with B ≤ L:")
    st.latex(r"s_c=1+0.3\frac{B}{L}")
    st.latex(r"s_q=1.0")
    st.latex(r"s_\gamma=1-0.2\frac{B}{L}")

with st.expander("4. Contact Pressure Under Axial Load and Biaxial Moment"):
    st.markdown("Coordinate convention: x is along footing width B, y is along footing length L.")
    st.latex(r"A=BL")
    st.latex(r"I_x=\frac{BL^3}{12}")
    st.latex(r"I_y=\frac{LB^3}{12}")
    st.latex(r"q(x,y)=\frac{P}{A}+\frac{M_yx}{I_y}+\frac{M_xy}{I_x}")
    st.markdown("Maximum and minimum pressures occur at footing corners for full-contact linear pressure distribution.")

with st.expander("5. Eccentricity and Kern Check"):
    st.latex(r"e_x=\frac{M_y}{P}")
    st.latex(r"e_y=\frac{M_x}{P}")
    st.latex(r"|e_x|\leq\frac{B}{6}")
    st.latex(r"|e_y|\leq\frac{L}{6}")
    st.markdown(
        "If eccentricity exceeds the kern limit, part of the footing may lose contact with the soil. "
        "A partial-contact analysis is then required."
    )

with st.expander("6. Groundwater Correction Used in This App"):
    st.markdown(
        "For the γBNγ term, the app applies a simplified groundwater correction:"
    )
    st.markdown(
        "- Water table below Df + B: use moist γ.\n"
        "- Water table at/above foundation base: use submerged γ = γsat − γw.\n"
        "- Water table between Df and Df + B: interpolate between moist and submerged γ."
    )
    st.latex(r"\gamma' = \gamma_{sat}-\gamma_w")

with st.expander("7. Important Design Limitations"):
    st.markdown(
        "This app checks bearing capacity only. A complete spread footing design should also evaluate:\n\n"
        "- Total and differential settlement\n"
        "- Sliding resistance\n"
        "- Overturning stability\n"
        "- Punching shear and one-way shear\n"
        "- Flexural reinforcement design\n"
        "- Minimum embedment, frost depth, expansive soil effects, scour, and drainage\n"
        "- Load combinations and strength/serviceability requirements from the applicable building code"
    )


# ============================================================
# Downloadable Results
# ============================================================
st.divider()
st.subheader("Export Calculation Summary")

export_df = pd.DataFrame(
    {
        "Item": [
            "P (kips)",
            "Mx (kip-ft)",
            "My (kip-ft)",
            "B (ft)",
            "L (ft)",
            "Df (ft)",
            "c (psf)",
            "phi (deg)",
            "gamma moist (pcf)",
            "gamma sat (pcf)",
            "water table depth (ft)",
            "FS target",
            "q avg gross (ksf)",
            "q max gross (ksf)",
            "q min gross (ksf)",
            "q surcharge (ksf)",
            "qult gross (ksf)",
            "qult net (ksf)",
            "qall gross (ksf)",
            "qall net (ksf)",
            "actual FS based on selected comparison",
            "status",
        ],
        "Value": [
            P_kips,
            Mx_kipft,
            My_kipft,
            B_ft,
            L_ft,
            Df_ft,
            c_psf,
            phi_deg,
            gamma_moist_pcf,
            gamma_sat_pcf,
            water_table_depth_ft,
            fs_target,
            pr.q_avg_ksf,
            pr.q_max_ksf,
            pr.q_min_ksf,
            bc.surcharge_q_ksf,
            bc.qult_gross_ksf,
            bc.qult_net_ksf,
            bc.qall_gross_ksf,
            bc.qall_net_ksf,
            fs_actual,
            status_title,
        ],
    }
)

csv_bytes = export_df.to_csv(index=False).encode("utf-8")
st.download_button(
    label="Download summary as CSV",
    data=csv_bytes,
    file_name="spread_footing_terzaghi_summary.csv",
    mime="text/csv",
)
