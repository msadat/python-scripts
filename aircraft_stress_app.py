import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from math import pi, sqrt


st.set_page_config(page_title="Aircraft Pavement Stress Calculator", layout="wide")

st.title("Aircraft Pavement Stress Calculator")
st.write(
    "Estimate vertical stress below pavement using individual aircraft wheel loads, "
    "finite tire contact area, and stress bulb superposition."
)


AIRCRAFT_DATA = {
    "Boeing 737-800": {
        "main_gear_load_lb": 130000,
        "tire_pressure_psi": 200,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },
    "Airbus A321": {
        "main_gear_load_lb": 150000,
        "tire_pressure_psi": 210,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },
    "Boeing 777-300ER": {
        "main_gear_load_lb": 550000,
        "tire_pressure_psi": 221,
        "wheel_coordinates_ft": [
            (-6.0, -10.0), (-6.0, -8.0),
            (0.0, -10.0), (0.0, -8.0),
            (6.0, -10.0), (6.0, -8.0),
            (-6.0, 8.0), (-6.0, 10.0),
            (0.0, 8.0), (0.0, 10.0),
            (6.0, 8.0), (6.0, 10.0),
        ],
    },
    "Airbus A380-800": {
        "main_gear_load_lb": 900000,
        "tire_pressure_psi": 220,
        "wheel_coordinates_ft": [
            (-6.0, -20.0), (-6.0, -18.0),
            (-2.0, -20.0), (-2.0, -18.0),
            (-6.0, 18.0), (-6.0, 20.0),
            (-2.0, 18.0), (-2.0, 20.0),
            (4.0, -5.0), (4.0, -3.0),
            (8.0, -5.0), (8.0, -3.0),
            (12.0, -5.0), (12.0, -3.0),
            (4.0, 3.0), (4.0, 5.0),
            (8.0, 3.0), (8.0, 5.0),
            (12.0, 3.0), (12.0, 5.0),
        ],
    },
}


def circular_loaded_area_stress_psf(q_psf, radius_ft, x_ft, y_ft, z_ft):
    """
    Approximate vertical stress from a uniformly loaded circular tire contact area.

    For directly beneath the center:
    sigma_z = q * [1 - 1 / (1 + (a/z)^2)^(3/2)]

    For off-center points, this uses an equivalent radial distance correction:
    z_eff = sqrt(z^2 + r^2)

    This is an educational approximation, not layered elastic pavement design.
    """

    r_ft = sqrt(x_ft**2 + y_ft**2)
    z_eff = sqrt(z_ft**2 + r_ft**2)

    if z_eff <= 0:
        return q_psf

    stress = q_psf * (
        1 - 1 / ((1 + (radius_ft / z_eff) ** 2) ** 1.5)
    )

    return stress


def calculate_wheel_stresses(
    wheel_load_lb,
    tire_pressure_psi,
    wheel_coordinates_ft,
    depths_ft,
    analysis_x_ft,
    analysis_y_ft,
):
    tire_pressure_psf = tire_pressure_psi * 144.0

    contact_area_sq_in = wheel_load_lb / tire_pressure_psi
    contact_area_sq_ft = contact_area_sq_in / 144.0
    contact_radius_ft = sqrt(contact_area_sq_ft / pi)

    rows = []

    for z in depths_ft:
        total_stress = 0.0

        for i, (x_wheel, y_wheel) in enumerate(wheel_coordinates_ft, start=1):
            x_offset = analysis_x_ft - x_wheel
            y_offset = analysis_y_ft - y_wheel

            stress_psf = circular_loaded_area_stress_psf(
                tire_pressure_psf,
                contact_radius_ft,
                x_offset,
                y_offset,
                z,
            )

            total_stress += stress_psf

            rows.append({
                "Depth Below Pavement (ft)": z,
                "Wheel": f"Wheel {i}",
                "Stress Contribution (psf)": stress_psf,
                "Stress Contribution (psi)": stress_psf / 144.0,
                "X Wheel Coordinate (ft)": x_wheel,
                "Y Wheel Coordinate (ft)": y_wheel,
            })

        rows.append({
            "Depth Below Pavement (ft)": z,
            "Wheel": "Total Combined Stress",
            "Stress Contribution (psf)": total_stress,
            "Stress Contribution (psi)": total_stress / 144.0,
            "X Wheel Coordinate (ft)": None,
            "Y Wheel Coordinate (ft)": None,
        })

    return pd.DataFrame(rows), contact_area_sq_in, contact_radius_ft


st.sidebar.header("Aircraft Selection")

aircraft = st.sidebar.selectbox("Select Aircraft", list(AIRCRAFT_DATA.keys()))
data = AIRCRAFT_DATA[aircraft]

wheel_coordinates = data["wheel_coordinates_ft"]
number_of_wheels = len(wheel_coordinates)

st.sidebar.header("Load Inputs")

total_main_gear_load = st.sidebar.number_input(
    "Total Modeled Main Gear Load (lb)",
    min_value=1000.0,
    value=float(data["main_gear_load_lb"]),
    step=1000.0,
)

equalized_wheel_load = total_main_gear_load / number_of_wheels

load_method = st.sidebar.radio(
    "Wheel Load Method",
    ["Equal load per tire", "User-defined tire load"],
)

if load_method == "Equal load per tire":
    wheel_load = equalized_wheel_load
else:
    wheel_load = st.sidebar.number_input(
        "Individual Tire Load (lb)",
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

st.sidebar.header("Analysis Settings")

max_depth = st.sidebar.number_input(
    "Maximum Depth (ft)",
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


depths = np.arange(depth_increment, max_depth + depth_increment, depth_increment)

stress_df, contact_area_sq_in, contact_radius_ft = calculate_wheel_stresses(
    wheel_load,
    tire_pressure,
    wheel_coordinates,
    depths,
    analysis_x,
    analysis_y,
)

total_df = stress_df[stress_df["Wheel"] == "Total Combined Stress"]
individual_df = stress_df[stress_df["Wheel"] != "Total Combined Stress"]

contact_radius_in = contact_radius_ft * 12.0
contact_diameter_in = contact_radius_in * 2.0


col1, col2, col3, col4 = st.columns(4)

col1.metric("Aircraft", aircraft)
col2.metric("Main Gear Wheels", number_of_wheels)
col3.metric("Individual Tire Load", f"{wheel_load:,.0f} lb")
col4.metric("Tire Pressure", f"{tire_pressure:,.0f} psi")


st.subheader("Tire Contact Area")

st.write(
    f"Each tire is modeled as a finite circular loaded area, not a point load.  \n"
    f"Approximate tire contact area = **{contact_area_sq_in:,.1f} in²**  \n"
    f"Equivalent contact radius = **{contact_radius_in:,.1f} in**  \n"
    f"Equivalent contact diameter = **{contact_diameter_in:,.1f} in**"
)


st.subheader("Combined Vertical Stress vs Depth")

fig_total = px.line(
    total_df,
    x="Depth Below Pavement (ft)",
    y="Stress Contribution (psf)",
    markers=True,
    title=f"Combined Stress Below Pavement - {aircraft}",
)

fig_total.update_layout(
    xaxis_title="Depth Below Pavement (ft)",
    yaxis_title="Combined Vertical Stress (psf)",
)

st.plotly_chart(fig_total, use_container_width=True)


st.subheader("Individual Wheel Stress Contributions and Overlap")

fig_overlap = px.line(
    stress_df,
    x="Depth Below Pavement (ft)",
    y="Stress Contribution (psf)",
    color="Wheel",
    title="Stress Contributions from Individual Wheels plus Total Combined Stress",
)

fig_overlap.update_layout(
    xaxis_title="Depth Below Pavement (ft)",
    yaxis_title="Vertical Stress (psf)",
)

st.plotly_chart(fig_overlap, use_container_width=True)


st.subheader("Stress Bulb Interaction at Selected Depth")

selected_depth = st.slider(
    "Select Depth for Plan-View Stress Bulb Map (ft)",
    min_value=float(depth_increment),
    max_value=float(max_depth),
    value=min(5.0, float(max_depth)),
    step=float(depth_increment),
)

grid_range = st.slider(
    "Plan View Grid Extent (ft)",
    min_value=10,
    max_value=80,
    value=40,
    step=5,
)

x_grid = np.linspace(-grid_range, grid_range, 80)
y_grid = np.linspace(-grid_range, grid_range, 80)

grid_rows = []

for x in x_grid:
    for y in y_grid:
        total_stress = 0.0

        for x_wheel, y_wheel in wheel_coordinates:
            x_offset = x - x_wheel
            y_offset = y - y_wheel

            stress_psf = circular_loaded_area_stress_psf(
                tire_pressure * 144.0,
                contact_radius_ft,
                x_offset,
                y_offset,
                selected_depth,
            )

            total_stress += stress_psf

        grid_rows.append({
            "X (ft)": x,
            "Y (ft)": y,
            "Stress (psf)": total_stress,
        })

grid_df = pd.DataFrame(grid_rows)

fig_contour = px.density_heatmap(
    grid_df,
    x="X (ft)",
    y="Y (ft)",
    z="Stress (psf)",
    nbinsx=80,
    nbinsy=80,
    title=f"Plan-View Combined Stress Bulb at {selected_depth:.1f} ft Depth",
)

wheel_df = pd.DataFrame(
    wheel_coordinates,
    columns=["X (ft)", "Y (ft)"]
)

fig_wheels = px.scatter(
    wheel_df,
    x="X (ft)",
    y="Y (ft)",
)

for trace in fig_wheels.data:
    fig_contour.add_trace(trace)

fig_contour.update_yaxes(scaleanchor="x", scaleratio=1)

st.plotly_chart(fig_contour, use_container_width=True)


st.subheader("Simplified Main Gear Layout")

fig_layout = px.scatter(
    wheel_df,
    x="X (ft)",
    y="Y (ft)",
    title=f"Simplified Main Gear Wheel Layout - {aircraft}",
)

fig_layout.update_traces(marker=dict(size=14))
fig_layout.update_yaxes(scaleanchor="x", scaleratio=1)

st.plotly_chart(fig_layout, use_container_width=True)


st.subheader("Combined Stress Table")

display_total_df = total_df[[
    "Depth Below Pavement (ft)",
    "Stress Contribution (psf)",
    "Stress Contribution (psi)",
]].rename(columns={
    "Stress Contribution (psf)": "Combined Vertical Stress (psf)",
    "Stress Contribution (psi)": "Combined Vertical Stress (psi)",
})

st.dataframe(
    display_total_df.style.format({
        "Depth Below Pavement (ft)": "{:.2f}",
        "Combined Vertical Stress (psf)": "{:,.1f}",
        "Combined Vertical Stress (psi)": "{:,.2f}",
    }),
    use_container_width=True,
)


st.warning(
    "Engineering note: This is an educational approximation. It models each aircraft tire "
    "as a finite circular loaded area and sums stress contributions from all wheels. "
    "It does not replace layered elastic analysis, FAARFIELD, gear wander analysis, "
    "pass-to-coverage calculation, or aircraft manufacturer pavement design data."
)
