import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from math import pi, sqrt


st.set_page_config(
    page_title="Aircraft Pavement Stress Calculator",
    layout="wide"
)

st.title("Aircraft Pavement Stress Calculator")
st.write(
    "Estimate vertical stress below pavement at different depths using "
    "simplified Boussinesq point-load superposition from individual aircraft wheel loads."
)


AIRCRAFT_DATA = {
    "Boeing 737-800": {
        "description": "Simplified dual-wheel main gear aircraft",
        "main_gear_load_lb": 130000,
        "tire_pressure_psi": 200,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },

    "Airbus A321": {
        "description": "Simplified dual-wheel main gear aircraft",
        "main_gear_load_lb": 150000,
        "tire_pressure_psi": 210,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },

    "Boeing 777-300ER": {
        "description": "Two 6-wheel main landing gear trucks, 12 main gear wheels total",
        "main_gear_load_lb": 550000,
        "tire_pressure_psi": 221,
        "wheel_coordinates_ft": [
            # Left main gear truck, 6 wheels
            (-6.0, -10.0), (-6.0, -8.0),
            (0.0, -10.0), (0.0, -8.0),
            (6.0, -10.0), (6.0, -8.0),

            # Right main gear truck, 6 wheels
            (-6.0, 8.0), (-6.0, 10.0),
            (0.0, 8.0), (0.0, 10.0),
            (6.0, 8.0), (6.0, 10.0),
        ],
    },

    "Airbus A380-800": {
        "description": "Simplified 20-wheel main landing gear configuration",
        "main_gear_load_lb": 900000,
        "tire_pressure_psi": 220,
        "wheel_coordinates_ft": [
            # Left wing gear, 4 wheels
            (-6.0, -20.0), (-6.0, -18.0),
            (-2.0, -20.0), (-2.0, -18.0),

            # Right wing gear, 4 wheels
            (-6.0, 18.0), (-6.0, 20.0),
            (-2.0, 18.0), (-2.0, 20.0),

            # Left body gear, 6 wheels
            (4.0, -5.0), (4.0, -3.0),
            (8.0, -5.0), (8.0, -3.0),
            (12.0, -5.0), (12.0, -3.0),

            # Right body gear, 6 wheels
            (4.0, 3.0), (4.0, 5.0),
            (8.0, 3.0), (8.0, 5.0),
            (12.0, 3.0), (12.0, 5.0),
        ],
    },
}


def boussinesq_vertical_stress_psf(Q_lb, x_ft, y_ft, z_ft):
    """
    Vertical stress beneath a point load using Boussinesq equation.

    sigma_z = 3Qz^3 / (2πR^5)

    Q = individual wheel load, lb
    x = horizontal x-offset from wheel, ft
    y = horizontal y-offset from wheel, ft
    z = depth below pavement surface, ft

    Output = vertical stress, psf
    """

    R = sqrt(x_ft**2 + y_ft**2 + z_ft**2)

    if R == 0:
        return 0.0

    sigma_z = (3 * Q_lb * z_ft**3) / (2 * pi * R**5)

    return sigma_z


def total_stress_at_depth(
    wheel_load_lb,
    wheel_coordinates_ft,
    z_ft,
    analysis_x_ft=0.0,
    analysis_y_ft=0.0,
):
    """
    Sum vertical stress contributions from each individual wheel.
    This is not a lumped aircraft load.
    """

    total_stress = 0.0

    for x_wheel, y_wheel in wheel_coordinates_ft:
        x_offset = analysis_x_ft - x_wheel
        y_offset = analysis_y_ft - y_wheel

        stress_from_wheel = boussinesq_vertical_stress_psf(
            wheel_load_lb,
            x_offset,
            y_offset,
            z_ft,
        )

        total_stress += stress_from_wheel

    return total_stress


st.sidebar.header("Aircraft Selection")

aircraft = st.sidebar.selectbox(
    "Select Aircraft",
    list(AIRCRAFT_DATA.keys())
)

data = AIRCRAFT_DATA[aircraft]
wheel_coordinates = data["wheel_coordinates_ft"]
number_of_wheels = len(wheel_coordinates)

st.sidebar.caption(data["description"])

st.sidebar.header("Load Assumptions")

total_main_gear_load = st.sidebar.number_input(
    "Total Modeled Main Gear Load (lb)",
    min_value=1000.0,
    value=float(data["main_gear_load_lb"]),
    step=1000.0,
)

equalized_wheel_load = total_main_gear_load / number_of_wheels

load_method = st.sidebar.radio(
    "Wheel Load Method",
    [
        "Equal load per main gear tire",
        "User-defined individual tire load",
    ],
)

if load_method == "Equal load per main gear tire":
    wheel_load = equalized_wheel_load
    st.sidebar.info(
        f"Wheel load = {total_main_gear_load:,.0f} lb / "
        f"{number_of_wheels} wheels = {wheel_load:,.0f} lb per tire"
    )
else:
    wheel_load = st.sidebar.number_input(
        "Individual Tire/Wheel Load (lb)",
        min_value=100.0,
        value=float(equalized_wheel_load),
        step=500.0,
    )

tire_pressure = st.sidebar.number_input(
    "Tire Pressure (psi)",
    min_value=50.0,
    value=float(data["tire_pressure_psi"]),
    step=5.0,
)

st.sidebar.header("Depth and Analysis Location")

max_depth = st.sidebar.number_input(
    "Maximum Depth Below Pavement (ft)",
    min_value=1.0,
    value=25.0,
    step=1.0,
)

depth_increment = st.sidebar.number_input(
    "Depth Increment (ft)",
    min_value=0.25,
    value=0.5,
    step=0.25,
)

analysis_x = st.sidebar.number_input(
    "Analysis Point X Coordinate (ft)",
    value=0.0,
    step=0.5,
)

analysis_y = st.sidebar.number_input(
    "Analysis Point Y Coordinate (ft)",
    value=0.0,
    step=0.5,
)


depths = np.arange(
    depth_increment,
    max_depth + depth_increment,
    depth_increment,
)

stress_values_psf = [
    total_stress_at_depth(
        wheel_load,
        wheel_coordinates,
        z,
        analysis_x,
        analysis_y,
    )
    for z in depths
]

stress_values_psi = [stress / 144.0 for stress in stress_values_psf]

df = pd.DataFrame({
    "Depth Below Pavement (ft)": depths,
    "Vertical Stress (psf)": stress_values_psf,
    "Vertical Stress (psi)": stress_values_psi,
})

contact_area_sq_in = wheel_load / tire_pressure
equivalent_contact_radius_in = sqrt(contact_area_sq_in / pi)
equivalent_contact_diameter_in = equivalent_contact_radius_in * 2


col1, col2, col3, col4 = st.columns(4)

col1.metric("Selected Aircraft", aircraft)
col2.metric("Main Gear Wheels Modeled", number_of_wheels)
col3.metric("Individual Wheel Load", f"{wheel_load:,.0f} lb")
col4.metric("Tire Pressure", f"{tire_pressure:,.0f} psi")


st.subheader("Load Distribution Method")

st.write(
    f"The model applies **{wheel_load:,.0f} lb** to each of the "
    f"**{number_of_wheels} individual main gear wheels**. "
    "Vertical stress is calculated separately for each wheel and then summed "
    "at the selected analysis point and depth."
)


st.subheader("Estimated Tire Contact Area")

st.write(
    f"Approximate single-tire contact area = **{contact_area_sq_in:,.1f} in²**  \n"
    f"Equivalent circular contact radius = **{equivalent_contact_radius_in:,.1f} in**  \n"
    f"Equivalent circular contact diameter = **{equivalent_contact_diameter_in:,.1f} in**"
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


st.subheader("Vertical Stress Table")

st.dataframe(
    df.style.format({
        "Depth Below Pavement (ft)": "{:.2f}",
        "Vertical Stress (psf)": "{:,.1f}",
        "Vertical Stress (psi)": "{:,.2f}",
    }),
    use_container_width=True,
)


st.subheader("Simplified Main Gear Wheel Layout")

layout_df = pd.DataFrame(
    wheel_coordinates,
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


st.subheader("Wheel Coordinate Table")

st.dataframe(
    layout_df,
    use_container_width=True,
)


st.warning(
    "Engineering note: This educational tool uses Boussinesq point-load "
    "superposition from individual wheel loads. It does not model layered elastic "
    "pavement response, tire contact pressure distribution, gear wander, "
    "load repetition, pass-to-coverage ratio, or FAARFIELD design methodology."
)
