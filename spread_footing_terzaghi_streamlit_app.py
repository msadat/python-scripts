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
# Data Classes
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


# ============================================================
# Bearing Capacity Functions
# ============================================================
def bearing_capacity_factors(phi_deg: float, ngamma_method: str) -> BearingCapacityFactors:
    """
    Classical bearing capacity factors.
    For phi = 0:
        Nc = 5.7, Nq = 1.0, Ngamma = 0.0
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
    """
    gamma_w = 62.4
    gamma_sub_pcf = max(gamma_sat_pcf - gamma_w, 0.0)

    if water_table_depth_ft <= foundation_depth_ft:
        return gamma_sub_pcf, "Water table at/above foundation base; submerged unit weight used for γ term."

    if water_table_depth_ft >= foundation_depth_ft + footing_width_ft:
        return gamma_moist_pcf, "Water table below Df + B; moist unit weight used for γ term."

    distance_below_base = water_table_depth_ft - foundation_depth_ft
    ratio_moist = max(0.0, min(distance_below_base / footing_width_ft, 1.0))
    gamma_eff = ratio_moist * gamma_moist_pcf + (1.0 - ratio_moist) * gamma_sub_pcf
    return gamma_eff, "Water table between Df and Df + B; interpolated effective γ used."


def surcharge_at_base(
    gamma_moist_pcf: float,
    gamma_sat_pcf: float,
    water_table_depth_ft: float,
    foundation_depth_ft: float,
) -> float:
    """
    Approximate effective surcharge at base, in ksf.
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
    Terzaghi-style rectangular footing shape factors.
    For B <= L:
        sc = 1 + 0.3(B/L)
        sq = 1.0
        sγ = 1 - 0.2(B/L)
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


# ============================================================
# Contact Pressure at Footing Base
# ============================================================
def compute_pressure_distribution(
    P_kips: float,
    Mx_kipft: float,
    My_kipft: float,
    B_ft: float,
    L_ft: float,
    grid_n: int = 81,
) -> PressureResult:
    """
    Rigid footing linear pressure distribution:
        q(x,y) = P/A + My*x/Iy + Mx*y/Ix

    x = width direction
    y = length direction

    Mx varies pressure along y.
    My varies pressure along x.
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


def compute_base_pressure_grid(
    P_kips: float,
    Mx_kipft: float,
    My_kipft: float,
    B_ft: float,
    L_ft: float,
    n: int = 41,
    clip_tension: bool = False,
):
    """
    Returns base contact pressure grid for use in subsurface stress calculation.
    """
    B = max(B_ft, 1e-9)
    L = max(L_ft, 1e-9)
    A = B * L
    Ix = B * L**3 / 12.0
    Iy = L * B**3 / 12.0

    x = np.linspace(-B / 2.0, B / 2.0, n)
    y = np.linspace(-L / 2.0, L / 2.0, n)
    X, Y = np.meshgrid(x, y)

    q_avg = P_kips / A if A > 0 else 0.0
    q = q_avg + (My_kipft * X / Iy if Iy > 0 else 0.0) + (Mx_kipft * Y / Ix if Ix > 0 else 0.0)

    if clip_tension:
        q = np.maximum(q, 0.0)

    return X, Y, q


# ============================================================
# Subsurface Stress (Boussinesq Numerical Integration)
# ============================================================
def compute_subsurface_stress_section(
    P_kips: float,
    Mx_kipft: float,
    My_kipft: float,
    B_ft: float,
    L_ft: float,
    z_max_ft: float,
    direction: str = "width",
    n_source: int = 31,
    n_eval: int = 51,
    lateral_extent_factor: float = 1.0,
    clip_tension: bool = True,
):
    """
    Computes vertical stress below the footing by numerically integrating the
    Boussinesq point-load solution over the loaded footing base.

    Vertical stress due to point load:
        σz = (3Q / 2π) * z^3 / R^5

    For distributed pressure:
        dQ = q dA
    """
    Xs, Ys, qs = compute_base_pressure_grid(
        P_kips, Mx_kipft, My_kipft, B_ft, L_ft, n=n_source, clip_tension=clip_tension
    )

    x_source = Xs[0, :]
    y_source = Ys[:, 0]
    dx = x_source[1] - x_source[0]
    dy = y_source[1] - y_source[0]

    if direction == "width":
        span = B_ft
        coord = np.linspace(-lateral_extent_factor * B_ft, lateral_extent_factor * B_ft, n_eval)
    else:
        span = L_ft
        coord = np.linspace(-lateral_extent_factor * L_ft, lateral_extent_factor * L_ft, n_eval)

    z_vals = np.linspace(0.1, z_max_ft, n_eval)
    sigma = np.zeros((len(z_vals), len(coord)))

    for iz, z in enumerate(z_vals):
        if direction == "width":
            # section at y = 0, varying x
            for ic, x0 in enumerate(coord):
                R2 = (x0 - Xs) ** 2 + Ys**2 + z**2
                sigma[iz, ic] = np.sum((3.0 * qs * z**3 / (2.0 * np.pi * (R2 ** 2.5))) * dx * dy)
        else:
            # section at x = 0, varying y
            for ic, y0 in enumerate(coord):
                R2 = Xs**2 + (y0 - Ys) ** 2 + z**2
                sigma[iz, ic] = np.sum((3.0 * qs * z**3 / (2.0 * np.pi * (R2 ** 2.5))) * dx * dy)

    return coord, z_vals, sigma


# ============================================================
# Plotting Functions
# ============================================================
def build_pressure_contour_figure(pr: PressureResult) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=pr.grid_x[0, :],
            y=pr.grid_y[:, 0],
            z=pr.grid_q,
            colorscale="Viridis",
            colorbar=dict(title="q (ksf)"),
            contours=dict(showlabels=True, labelfont=dict(size=10)),
        )
    )
    fig.update_layout(
        title="Plan View Contact Pressure at Footing Base",
        xaxis_title="Width direction, x (ft)",
        yaxis_title="Length direction, y (ft)",
        height=520,
    )
    return fig


def build_pressure_surface_figure(pr: PressureResult) -> go.Figure:
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
        title="3D Contact Pressure Surface",
        scene=dict(
            xaxis_title="x (ft)",
            yaxis_title="y (ft)",
            zaxis_title="q (ksf)",
        ),
        margin=dict(l=0, r=0, b=0, t=45),
        height=600,
    )
    return fig


def build_base_profile_figure(coord: np.ndarray, qline: np.ndarray, axis_name: str, moment_name: str) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=coord,
            y=qline,
            mode="lines",
            fill="tozeroy",
            name="Base pressure",
        )
    )
    fig.update_layout(
        title=f"2D Base Pressure Profile Along {axis_name} Direction ({moment_name} effect)",
        xaxis_title=f"{axis_name} coordinate (ft)",
        yaxis_title="Bearing pressure at base, q (ksf)",
        height=420,
    )
    return fig


def build_subsurface_heatmap_figure(
    coord: np.ndarray,
    z_vals: np.ndarray,
    sigma: np.ndarray,
    section_name: str,
) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            x=coord,
            y=z_vals,
            z=sigma,
            colorscale="Viridis",
            colorbar=dict(title="σz (ksf)"),
        )
    )
    fig.update_layout(
        title=f"Subsurface Vertical Stress Distribution - {section_name} Section",
        xaxis_title=f"{section_name} coordinate (ft)",
        yaxis_title="Depth below footing base, z (ft)",
        height=500,
        yaxis=dict(autorange="reversed"),
    )
    return fig


def build_centerline_depth_figure(z_vals: np.ndarray, sigma_centerline: np.ndarray, section_name: str) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=sigma_centerline,
            y=z_vals,
            mode="lines",
            name="σz",
        )
    )
    fig.update_layout(
        title=f"Vertical Stress vs Depth at Centerline - {section_name} Section",
        xaxis_title="Vertical stress, σz (ksf)",
        yaxis_title="Depth below footing base, z (ft)",
        height=420,
        yaxis=dict(autorange="reversed"),
    )
    return fig


def build_loading_schematic(
    span_ft: float,
    Df_ft: float,
    coord_line: np.ndarray,
    q_line: np.ndarray,
    section_title: str,
    section_axis: str,
    moment_label: str,
) -> go.Figure:
    """
    Simple 2D schematic showing footing, axial load, applied moment,
    and pressure profile under the footing.
    """
    thickness = max(0.8, 0.12 * span_ft)
    col_w = 0.18 * span_ft
    col_h = 0.35 * span_ft

    y_ground = 0.0
    y_base_top = -Df_ft
    y_base_bot = -Df_ft - thickness

    q_comp = np.maximum(q_line, 0.0)
    qmax = max(np.max(q_comp), 1e-6)
    pressure_scale = 0.35 * span_ft / qmax
    y_pressure = y_base_bot - q_comp * pressure_scale

    fig = go.Figure()

    # Ground line
    fig.add_shape(type="line", x0=-0.9 * span_ft, x1=0.9 * span_ft, y0=y_ground, y1=y_ground)

    # Footing
    fig.add_shape(
        type="rect",
        x0=-span_ft / 2.0,
        x1=span_ft / 2.0,
        y0=y_base_bot,
        y1=y_base_top,
        line=dict(color="black"),
        fillcolor="lightgray",
    )

    # Column / pedestal
    fig.add_shape(
        type="rect",
        x0=-col_w / 2.0,
        x1=col_w / 2.0,
        y0=y_base_top,
        y1=y_base_top + col_h,
        line=dict(color="black"),
        fillcolor="silver",
    )

    # Axial load arrow
    fig.add_annotation(
        x=0.0,
        y=y_base_top + col_h,
        ax=0.0,
        ay=y_base_top + col_h + 0.35 * span_ft,
        showarrow=True,
        arrowhead=3,
        arrowsize=1.4,
        arrowwidth=2,
        text="P",
    )

    # Moment arc
    theta = np.linspace(0.2 * np.pi, 1.45 * np.pi, 60)
    r = 0.22 * span_ft
    x_arc = r * np.cos(theta)
    y_arc = (y_base_top + col_h + 0.12 * span_ft) + r * np.sin(theta)
    fig.add_trace(go.Scatter(x=x_arc, y=y_arc, mode="lines", name=moment_label))
    fig.add_annotation(
        x=x_arc[-1],
        y=y_arc[-1],
        ax=x_arc[-5],
        ay=y_arc[-5],
        showarrow=True,
        arrowhead=2,
        arrowsize=1.2,
        arrowwidth=2,
        text=moment_label,
    )

    # Pressure polygon
    poly_x = list(coord_line) + list(coord_line[::-1])
    poly_y = list(y_pressure) + [y_base_bot] * len(coord_line)
    fig.add_trace(
        go.Scatter(
            x=poly_x,
            y=poly_y,
            fill="toself",
            mode="lines",
            line=dict(color="blue"),
            fillcolor="rgba(0, 0, 255, 0.25)",
            name="Bearing pressure",
        )
    )

    fig.add_annotation(
        x=coord_line[int(0.8 * len(coord_line))],
        y=np.min(y_pressure),
        text="q",
        showarrow=False,
        font=dict(color="blue"),
    )

    fig.add_annotation(
        x=-0.82 * span_ft,
        y=y_ground + 0.05 * span_ft,
        text="Ground surface",
        showarrow=False,
    )

    fig.add_annotation(
        x=0.65 * span_ft,
        y=(y_base_top + y_base_bot) / 2.0,
        text=f"Df = {Df_ft:.2f} ft",
        showarrow=False,
    )

    fig.update_layout(
        title=section_title,
        xaxis_title=f"{section_axis} direction (ft)",
        yaxis_title="Elevation / schematic depth",
        height=520,
        showlegend=False,
    )

    fig.update_yaxes(scaleanchor=None)
    return fig


# ============================================================
# Status Helper
# ============================================================
def classify_status(fs_actual: float, fs_target: float, uplift: bool) -> Tuple[str, str]:
    if uplift:
        return (
            "⚠️ Not adequate for full-contact assumption",
            "Negative bearing pressure is predicted over part of the footing. Soil cannot resist tension. Increase footing size, reduce eccentricity, or use partial-contact/effective-area analysis."
        )

    if fs_actual >= fs_target:
        return (
            "✅ Adequate for bearing pressure",
            "The footing satisfies the selected bearing capacity check."
        )

    return (
        "❌ Not adequate for bearing capacity",
        "The footing does not satisfy the selected bearing pressure check. Increase footing dimensions, improve soil, increase embedment, or reduce loads."
    )


# ============================================================
# Header
# ============================================================
st.title("Spread Footing Design Using Terzaghi Bearing Capacity Theory")
st.caption(
    "Preliminary geotechnical design tool for rectangular spread footings under axial load and biaxial moments. "
    "Terzaghi is used for bearing capacity; Boussinesq elastic integration is used to visualize subsurface stress spread."
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
    "Depth to groundwater below ground surface (ft)",
    min_value=0.0,
    value=999.0,
    step=1.0,
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
grid_n = st.sidebar.slider("Base pressure plot resolution", min_value=25, max_value=151, value=81, step=2)
z_max_ft = st.sidebar.number_input("Maximum depth below footing base for subsurface plots (ft)", min_value=1.0, value=15.0, step=1.0)
n_eval = st.sidebar.slider("Subsurface plot resolution", min_value=21, max_value=81, value=41, step=2)
n_source = st.sidebar.slider("Subsurface integration grid", min_value=15, max_value=51, value=31, step=2)
lateral_extent_factor = st.sidebar.slider("Subsurface lateral extent factor", min_value=0.5, max_value=2.0, value=1.0, step=0.1)


# ============================================================
# Calculations
# ============================================================
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
    gamma_moist_pcf,
    gamma_sat_pcf,
    water_table_depth_ft,
    Df_ft,
    min(B_ft, L_ft),
)

service_net_qmax = max(pr.q_max_ksf - bc.surcharge_q_ksf, 0.0)

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
status_title, status_body = classify_status(fs_actual, fs_target, pr.uplift_area_present)

# Base pressure profile lines
x_line = pr.grid_x[0, :]
y_line = pr.grid_y[:, 0]
row_mid = len(y_line) // 2
col_mid = len(x_line) // 2

q_width_line = pr.grid_q[row_mid, :]   # at y = 0
q_length_line = pr.grid_q[:, col_mid]  # at x = 0

# Subsurface stress sections using Boussinesq integration
coord_w, z_vals_w, sigma_w = compute_subsurface_stress_section(
    P_kips=P_kips,
    Mx_kipft=Mx_kipft,
    My_kipft=My_kipft,
    B_ft=B_ft,
    L_ft=L_ft,
    z_max_ft=z_max_ft,
    direction="width",
    n_source=n_source,
    n_eval=n_eval,
    lateral_extent_factor=lateral_extent_factor,
    clip_tension=True,
)

coord_l, z_vals_l, sigma_l = compute_subsurface_stress_section(
    P_kips=P_kips,
    Mx_kipft=Mx_kipft,
    My_kipft=My_kipft,
    B_ft=B_ft,
    L_ft=L_ft,
    z_max_ft=z_max_ft,
    direction="length",
    n_source=n_source,
    n_eval=n_eval,
    lateral_extent_factor=lateral_extent_factor,
    clip_tension=True,
)

sigma_w_center = sigma_w[:, len(coord_w) // 2]
sigma_l_center = sigma_l[:, len(coord_l) // 2]


# ============================================================
# Main Summary
# ============================================================
st.subheader("Design Summary")

c1, c2, c3, c4 = st.columns(4)
c1.metric("Footing Area", f"{A_ft2:,.2f} ft²")
c2.metric("Average Base Pressure", f"{pr.q_avg_ksf:,.3f} ksf")
c3.metric("Maximum Base Pressure", f"{pr.q_max_ksf:,.3f} ksf")
c4.metric("Minimum Base Pressure", f"{pr.q_min_ksf:,.3f} ksf")

box = st.container(border=True)
with box:
    st.markdown(f"### {status_title}")
    st.write(status_body)
    st.write(
        f"**Selected comparison basis:** {basis_label.capitalize()} pressure"
        f"  |  Demand = `{demand_ksf:,.3f} ksf`"
        f"  |  Allowable = `{allowable_ksf:,.3f} ksf`"
        f"  |  FS = `{fs_actual:,.3f}`"
        f"  |  Target FS = `{fs_target:,.3f}`"
    )

if not pr.kern_ok_x or not pr.kern_ok_y:
    st.warning(
        "The resultant is outside the kern in at least one direction. "
        "This indicates possible loss of contact over part of the footing."
    )

if pr.uplift_area_present:
    st.error(
        "Negative contact pressure exists in the full-contact solution. "
        "For subsurface stress plots, tensile base pressure has been clipped to zero for visualization."
    )


# ============================================================
# Detailed Tables
# ============================================================
st.subheader("Bearing Capacity Calculation")

calc_df = pd.DataFrame({
    "Parameter": [
        "Nc",
        "Nq",
        "Nγ",
        "Shape factor sc",
        "Shape factor sq",
        "Shape factor sγ",
        "Surcharge at base, q",
        "Effective γ for bearing term",
        "Ultimate gross bearing capacity, qult,gross",
        "Ultimate net bearing capacity, qult,net",
        "Allowable gross bearing pressure, qall,gross",
        "Allowable net bearing pressure, qall,net",
        "e_x = My / P",
        "e_y = Mx / P",
        "B / 6",
        "L / 6",
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
        f"{pr.e_x_ft:,.3f} ft",
        f"{pr.e_y_ft:,.3f} ft",
        f"{B_ft/6.0:,.3f} ft",
        f"{L_ft/6.0:,.3f} ft",
    ],
})
st.dataframe(calc_df, use_container_width=True, hide_index=True)
st.info(gw_note)


# ============================================================
# Tabs for Graphics
# ============================================================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Plan & Base Profiles",
    "Subsurface Stress",
    "Foundation Schematics",
    "3D Pressure Surface",
    "Formula Reference",
])

with tab1:
    st.subheader("2D Contact Pressure Views")

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(build_pressure_contour_figure(pr), use_container_width=True)
    with col2:
        st.write(
            "This contour shows the **contact pressure at the footing base**. "
            "Pressure varies linearly because the footing is treated as rigid under axial load and biaxial moment."
        )
        st.plotly_chart(build_base_profile_figure(x_line, q_width_line, "Width", "My"), use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(build_base_profile_figure(y_line, q_length_line, "Length", "Mx"), use_container_width=True)
    with col4:
        st.write(
            "The width profile is taken at **y = 0**, and the length profile is taken at **x = 0**. "
            "These are useful 2D section views of the pressure distribution under the footing."
        )

with tab2:
    st.subheader("Subsurface Stress Variation with Depth and Width / Length")

    st.markdown(
        """
        The plots below show **vertical stress in the soil mass below the footing**.  
        These are computed using **numerical integration of the Boussinesq point-load solution**
        over the base contact pressure distribution.

        **Important:**  
        - This is **not** the Terzaghi bearing capacity equation.  
        - Terzaghi is used for **ultimate bearing capacity**.  
        - Boussinesq here is used to **visualize how stress spreads with depth**.
        """
    )

    r1c1, r1c2 = st.columns(2)
    with r1c1:
        st.plotly_chart(
            build_subsurface_heatmap_figure(coord_w, z_vals_w, sigma_w, "Width"),
            use_container_width=True,
        )
    with r1c2:
        st.plotly_chart(
            build_centerline_depth_figure(z_vals_w, sigma_w_center, "Width"),
            use_container_width=True,
        )

    r2c1, r2c2 = st.columns(2)
    with r2c1:
        st.plotly_chart(
            build_subsurface_heatmap_figure(coord_l, z_vals_l, sigma_l, "Length"),
            use_container_width=True,
        )
    with r2c2:
        st.plotly_chart(
            build_centerline_depth_figure(z_vals_l, sigma_l_center, "Length"),
            use_container_width=True,
        )

with tab3:
    st.subheader("2D Foundation Pictures Showing Axial Load and Moments")

    s1, s2 = st.columns(2)
    with s1:
        st.plotly_chart(
            build_loading_schematic(
                span_ft=B_ft,
                Df_ft=Df_ft,
                coord_line=x_line,
                q_line=q_width_line,
                section_title="Width Section: Axial Load P and Moment My",
                section_axis="Width, x",
                moment_label="My",
            ),
            use_container_width=True,
        )

    with s2:
        st.plotly_chart(
            build_loading_schematic(
                span_ft=L_ft,
                Df_ft=Df_ft,
                coord_line=y_line,
                q_line=q_length_line,
                section_title="Length Section: Axial Load P and Moment Mx",
                section_axis="Length, y",
                moment_label="Mx",
            ),
            use_container_width=True,
        )

    st.write(
        "These schematics are conceptual 2D section views. "
        "They graphically show the applied axial load, the corresponding moment direction, "
        "and the resulting pressure shape beneath the footing."
    )

with tab4:
    st.subheader("3D Contact Pressure Surface")
    st.plotly_chart(build_pressure_surface_figure(pr), use_container_width=True)

with tab5:
    st.header("Formula Reference")

    with st.expander("1. Terzaghi Bearing Capacity Equation", expanded=True):
        st.latex(r"q_{ult,gross}=cN_c s_c + qN_q s_q + \frac{1}{2}\gamma B N_\gamma s_\gamma")
        st.markdown(
            "`B` is the smaller footing dimension in the bearing term, "
            "`q` is the surcharge at foundation base, and `sc`, `sq`, `sγ` are shape factors."
        )
        st.latex(r"q_{ult,net}=q_{ult,gross}-q")
        st.latex(r"q_{all}=\frac{q_{ult}}{FS}")

    with st.expander("2. Bearing Capacity Factors"):
        st.latex(r"N_q=e^{\pi \tan\phi}\tan^2\left(45^\circ+\frac{\phi}{2}\right)")
        st.latex(r"N_c=\frac{N_q-1}{\tan\phi}")
        st.markdown("For φ = 0 undrained clay:")
        st.latex(r"N_c=5.7,\quad N_q=1.0,\quad N_\gamma=0.0")
        st.markdown("Selectable Nγ equations:")
        st.latex(r"N_\gamma=1.5(N_q-1)\tan\phi \quad \text{Terzaghi approximation}")
        st.latex(r"N_\gamma=(N_q-1)\tan(1.4\phi) \quad \text{Meyerhof}")
        st.latex(r"N_\gamma=2(N_q+1)\tan\phi \quad \text{Vesic}")

    with st.expander("3. Shape Factors for Rectangular Footing"):
        st.latex(r"s_c=1+0.3\frac{B}{L}")
        st.latex(r"s_q=1.0")
        st.latex(r"s_\gamma=1-0.2\frac{B}{L}")

    with st.expander("4. Contact Pressure Under Axial Load and Biaxial Moment"):
        st.latex(r"A=BL")
        st.latex(r"I_x=\frac{BL^3}{12}")
        st.latex(r"I_y=\frac{LB^3}{12}")
        st.latex(r"q(x,y)=\frac{P}{A}+\frac{M_yx}{I_y}+\frac{M_xy}{I_x}")
        st.markdown(
            "Pressure varies linearly across the footing for a rigid footing with full contact."
        )

    with st.expander("5. Eccentricity / Kern Check"):
        st.latex(r"e_x=\frac{M_y}{P}")
        st.latex(r"e_y=\frac{M_x}{P}")
        st.latex(r"|e_x|\leq\frac{B}{6}")
        st.latex(r"|e_y|\leq\frac{L}{6}")
        st.markdown(
            "If eccentricity exceeds the kern limit, part of the footing may theoretically lift off and a full-contact pressure assumption is no longer valid."
        )

    with st.expander("6. Boussinesq Vertical Stress Used for Subsurface Plots"):
        st.markdown("For a point load Q applied at the ground/foundation interface:")
        st.latex(r"\sigma_z=\frac{3Q}{2\pi}\frac{z^3}{R^5}")
        st.markdown("For distributed pressure, the app numerically integrates over the footing base:")
        st.latex(r"dQ=q\,dA")
        st.latex(r"\sigma_z=\iint \frac{3\,q\,z^3}{2\pi R^5}\,dA")
        st.markdown(
            "This is used only for the subsurface stress visualization, not for ultimate bearing capacity."
        )

    with st.expander("7. Limitations"):
        st.markdown(
            """
            This app is for preliminary design and educational use. A complete foundation design should also include:
            - total and differential settlement
            - sliding
            - overturning
            - structural flexure
            - one-way shear and punching shear
            - code-required load combinations
            - frost, expansive soil, drainage, and construction considerations
            """
        )


# ============================================================
# Export
# ============================================================
st.divider()
st.subheader("Export Summary")

export_df = pd.DataFrame({
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
        "q_avg base (ksf)",
        "q_max base (ksf)",
        "q_min base (ksf)",
        "surcharge q (ksf)",
        "qult gross (ksf)",
        "qult net (ksf)",
        "qall gross (ksf)",
        "qall net (ksf)",
        "actual FS",
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
    ]
})

csv_bytes = export_df.to_csv(index=False).encode("utf-8")
st.download_button(
    label="Download summary as CSV",
    data=csv_bytes,
    file_name="spread_footing_terzaghi_summary.csv",
    mime="text/csv",
)
