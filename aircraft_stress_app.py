import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from math import pi, sqrt

st.set_page_config(page_title="Aircraft Pavement Stress Calculator", layout="wide")

st.title("Aircraft Pavement Stress Calculator")
st.write(
    "Estimate vertical stress below pavement from simplified aircraft wheel loads "
    "using Boussinesq point-load superposition."
)

# -----------------------------
# Default aircraft assumptions
# -----------------------------
AIRCRAFT_DATA = {
    "Boeing 737-800": {
        "main_gear_load_lb": 130000,
        "wheel_load_lb": 32500,
        "wheel_coordinates_ft": [(-1.5, -1.0), (-1.5, 1.0), (1.5, -1.0), (1.5, 1.0)],
        "tire_pressure_psi": 200,
    },
    "Airbus A321": {
        "main_gear_load_lb": 150000,
        "wheel_load_lb": 37500,
        "wheel_coordinates_ft": [(-1.5, -1.0), (-1.5, 1.0), (1.5, -1.0), (1.5, 1.0)],
        "tire_pressure_psi": 210,
    },
    "Boeing 777-300ER": {
        "main_gear_load_lb": 550000,
        "wheel_load_lb": 45833,
        "wheel_coordinates_ft": [
            (-6, -2), (-6, 0), (-6, 2),
            (0, -2), (0, 0), (0, 2),
            (6, -2), (6, 0), (6, 2),
            (12, -2), (12, 0), (12, 2),
        ],
        "tire_pressure_psi": 221,
    },
    "Airbus A380-800": {
        "main_gear_load_lb": 850000,
        "wheel_load_lb": 35417,
        "wheel_coordinates_ft": [
            (-10, -6), (-10, -4), (-10, 4), (-10, 6),
            (-4, -6), (-4, -4), (-4, 4), (-4, 6),
            (4, -6), (4, -4), (4, 4), (4, 6),
            (10, -6), (10, -4), (10, 4), (10, 6),
            (16, -6), (16, -4), (16, 4), (16, 6),
            (22, -6), (22, -4), (22, 4), (22, 6),
        ],
        "tire_pressure_psi": 220,
    },
}


def boussinesq_vertical_stress_psf(Q_lb, x_ft, y_ft, z_ft):
    """
    Vertical stress below a point load using Boussinesq equation.

    sigma_z = 3Qz^3 / (2πR^5)

    Q = wheel load, lb
    x, y = horizontal offset from wheel, ft
    z = depth, ft
    output = psf
    """
    R = sqrt(x_ft**2 + y_ft**2 + z_ft**2)
    if R == 0:
        return 0
    return (3 * Q_lb * z_ft**3) / (2 * pi * R**5)


def total_stress_at_depth(wheel_load_lb, wheel_coordinates, z_ft, analysis_x=0, analysis_y=0):
    total = 0
    for xw, yw in wheel_coordinates:
        x_offset = analysis_x - xw
        y_offset = analysis_y - yw
        total += boussinesq_vertical_stress_psf(wheel_load_lb, x_offset, y_offset, z_ft)
    return total


# -----------------------------
# Sidebar inputs
# -----------------------------
st.sidebar.header("Aircraft Selection")

aircraft = st.sidebar.selectbox("Select Aircraft", list(AIRCRAFT_DATA.keys()))
data = AIRCRAFT_DATA[aircraft]

st.sidebar.header("Editable Load Assumptions")

main_gear_load = st.sidebar.number_input(
    "Total Main Gear Load (lb)",
    min_value=1000.0,
    value=float(data["main_gear_load_lb"]),
    step=1000.0,
)

number_of_wheels = len(data["wheel_coordinates_ft"])

wheel_load = st.sidebar.number_input(
    "Single Wheel Load (lb)",
    min_value=100.0,
    value=float(main_gear_load / number_of_wheels),
    step=500.0,
)

tire_pressure = st.sidebar.number_input(
    "Tire Pressure (psi)",
    min_value=50.0,
    value=float(data["tire_pressure_psi"]),
    step=5.0,
)

max_depth = st.sidebar.number_input("Maximum Depth (ft)", min_value=1.0, value=25.0, step=1.0)
depth_increment = st.sidebar.number_input("Depth Increment (ft)", min_value=0.25, value=0.5, step=0.25)

analysis_x = st.sidebar.number_input("Analysis Point X Offset (ft)", value=0.0, step=0.5)
analysis_y = st.sidebar.number_input("Analysis Point Y Offset (ft)", value=0.0, step=0.5)

# -----------------------------
# Calculations
# -----------------------------
depths = np.arange(depth_increment, max_depth + depth_increment, depth_increment)

stress_values = [
    total_stress_at_depth(
        wheel_load,
        data["wheel_coordinates_ft"],
        z,
        analysis_x,
        analysis_y,
    )
    for z in depths
]

stress_psi = [s / 144 for s in stress_values]

df = pd.DataFrame({
    "Depth Below Pavement (ft)": depths,
    "Vertical Stress (psf)": stress_values,
    "Vertical Stress (psi)": stress_psi,
})

contact_area_sq_in = wheel_load / tire_pressure
equivalent_contact_radius_in = sqrt(contact_area_sq_in / pi)

# -----------------------------
# Output
# -----------------------------
col1, col2, col3, col4 = st.columns(4)

col1.metric("Aircraft", aircraft)
col2.metric("No. of Main Gear Wheels Modeled", number_of_wheels)
col3.metric("Wheel Load", f"{wheel_load:,.0f} lb")
col4.metric("Tire Pressure", f"{tire_pressure:,.0f} psi")

st.subheader("Estimated Tire Contact Area")
st.write(
    f"Approximate single tire contact area = **{contact_area_sq_in:,.1f} in²**; "
    f"equivalent circular contact radius = **{equivalent_contact_radius_in:,.1f} in**."
)

st.subheader("Vertical Stress vs Depth")

fig = px.line(
    df,
    x="Depth Below Pavement (ft)",
    y="Vertical Stress (psf)",
    markers=True,
    title=f"Estimated Vertical Stress Below Pavement - {aircraft}",
)

fig.update_layout(
    xaxis_title="Depth Below Pavement (ft)",
    yaxis_title="Vertical Stress (psf)",
)

st.plotly_chart(fig, use_container_width=True)

st.subheader("Stress Table")
st.dataframe(df, use_container_width=True)

st.subheader("Simplified Wheel Layout")
layout_df = pd.DataFrame(
    data["wheel_coordinates_ft"],
    columns=["X Coordinate (ft)", "Y Coordinate (ft)"]
)

fig_layout = px.scatter(
    layout_df,
    x="X Coordinate (ft)",
    y="Y Coordinate (ft)",
    title=f"Simplified Main Gear Wheel Layout - {aircraft}",
)

fig_layout.update_traces(marker=dict(size=14))
fig_layout.update_yaxes(scaleanchor="x", scaleratio=1)

st.plotly_chart(fig_layout, use_container_width=True)

st.warning(
    "Engineering note: This tool uses simplified Boussinesq point-load superposition. "
    "It does not account for layered pavement stiffness, tire contact pressure distribution, "
    "wander, pass-to-coverage ratio, gear interaction in layered media, nonlinear subgrade behavior, "
    "or FAA FAARFIELD design methodology. Use for educational screening only."
)
