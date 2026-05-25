
import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from math import pi, sqrt


st.set_page_config(page_title="Aircraft Pavement Stress Calculator", layout="wide")

st.title("Aircraft Pavement Stress Calculator")
st.write(
    "Educational pavement stress calculator using individual aircraft wheel loads, "
    "fixed tire contact area, calculated contact pressure, and stress bulb superposition."
)


AIRCRAFT_DATA = {
    "Boeing 737-800": {
        "main_gear_load_lb": 130000,
        "contact_radius_ft": 0.65,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },
    "Airbus A321": {
        "main_gear_load_lb": 150000,
        "contact_radius_ft": 0.70,
        "wheel_coordinates_ft": [
            (-3.0, -5.0), (-3.0, -3.0),
            (-3.0, 3.0), (-3.0, 5.0),
        ],
    },
    "Boeing 777-300ER": {
        "main_gear_load_lb": 550000,
        "contact_radius_ft": 0.85,
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
        "contact_radius_ft": 0.90,
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

    Directly beneath center:
    sigma_z = q * [1 - 1 / (1 + (a/z)^2)^(3/2)]

    Off-center approximation:
    z_eff = sqrt(z^2 + r^2)

    q_psf = calculated tire-pavement contact pressure
    radius_ft = fixed equivalent tire contact radius
    """

    r_ft = sqrt(x_ft**2 + y_ft**2)
    z_eff = sqrt(z_ft**2 + r_ft**2)

    if z_eff <= 0:
        return q_psf

    return q_psf * (1 - 1 / ((1 + (radius_ft / z_eff) ** 2) ** 1.5))


def gear_centroid(wheel_coordinates):
    x_vals = [x for x, y in wheel_coordinates]
    y_vals = [y for x, y in wheel_coordinates]
    return sum(x_vals) / len(x_vals), sum(y_vals) / len(y_vals)


def total_stress_at_point(
    wheel_coordinates,
    contact_pressure_psf,
    contact_radius_ft,
    analysis_x,
    analysis_y,
    depth,
):
    total = 0.0

    for xw, yw in wheel_coordinates:
        total += circular_loaded_area_stress_psf(
            contact_pressure_psf,
            contact_radius_ft,
            analysis_x - xw,
            analysis_y - yw,
            depth,
        )

    return total


def single_wheel_center_stress(contact_pressure_psf, contact_radius_ft, depth):
    return circular_loaded_area_stress_psf(
        contact_pressure_psf,
        contact_radius_ft,
        0.0,
        0.0,
        depth,
    )


st.sidebar.header("Aircraft Selection")

aircraft = st.sidebar.selectbox("Select Aircraft", list(AIRCRAFT_DATA.keys()))
data = AIRCRAFT_DATA[aircraft]

wheel_coordinates = data["wheel_coordinates_ft"]
number_of_wheels = len(wheel_coordinates)

centroid_x, centroid_y = gear_centroid(wheel_coordinates)

st.sidebar.header("Load and Tire Inputs")

total_main_gear_load = st.sidebar.number_input(
    "Total Modeled Main Gear Load (lb)",
    min_value=1000.0,
    value=float(data["main_gear_load_lb"]),
    step=1000.0,
)

wheel_load = total_main_gear_load / number_of_wheels

contact_radius_ft = st.sidebar.number_input(
    "Fixed Equivalent Tire Contact Radius (ft)",
    min_value=0.10,
    value=float(data["contact_radius_ft"]),
    step=0.05,
)

contact_area_sq_ft = pi * contact_radius_ft**2
contact_area_sq_in = contact_area_sq_ft * 144.0

contact_pressure_psf = wheel_load / contact_area_sq_ft
contact_pressure_psi = contact_pressure_psf / 144.0

st.sidebar.caption(
    f"Individual tire load = {total_main_gear_load:,.0f} lb / "
    f"{number_of_wheels} wheels = {wheel_load:,.0f} lb per tire."
)

st.sidebar.caption(
    f"Calculated contact pressure = {wheel_load:,.0f} lb / "
    f"{contact_area_sq_in:,.1f} in² = {contact_pressure_psi:,.1f} psi."
)

st.sidebar.header("Depth Settings")

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

st.sidebar.header("Analysis Location")

analysis_mode = st.sidebar.radio(
    "Select Analysis Location",
    [
        "Under first wheel center",
        "Under gear centroid",
        "User-defined point",
    ],
)

if analysis_mode == "Under first wheel center":
    analysis_x, analysis_y = wheel_coordinates[0]
elif analysis_mode == "Under gear centroid":
    analysis_x, analysis_y = centroid_x, centroid_y
else:
    analysis_x = st.sidebar.number_input("Analysis X Coordinate (ft)", value=0.0, step=0.5)
    analysis_y = st.sidebar.number_input("Analysis Y Coordinate (ft)", value=0.0, step=0.5)


depths = np.arange(depth_increment, max_depth + depth_increment, depth_increment)

stress_rows = []

for z in depths:
    single = single_wheel_center_stress(contact_pressure_psf, contact_radius_ft, z)

    total = total_stress_at_point(
        wheel_coordinates,
        contact_pressure_psf,
        contact_radius_ft,
        analysis_x,
        analysis_y,
        z,
    )

    stress_rows.append({
        "Depth Below Pavement (ft)": z,
        "Stress Type": "Single Wheel Centerline",
        "Vertical Stress (psf)": single,
        "Vertical Stress (psi)": single / 144.0,
    })

    stress_rows.append({
        "Depth Below Pavement (ft)": z,
        "Stress Type": "Combined Gear at Selected Point",
        "Vertical Stress (psf)": total,
        "Vertical Stress (psi)": total / 144.0,
    })

comparison_df = pd.DataFrame(stress_rows)

wheel_contribution_rows = []

for z in depths:
    for i, (xw, yw) in enumerate(wheel_coordinates, start=1):
        stress = circular_loaded_area_stress_psf(
            contact_pressure_psf,
            contact_radius_ft,
            analysis_x - xw,
            analysis_y - yw,
            z,
        )

        wheel_contribution_rows.append({
            "Depth Below Pavement (ft)": z,
            "Wheel": f"Wheel {i}",
            "Stress Contribution (psf)": stress,
            "Stress Contribution (psi)": stress / 144.0,
            "Distance From Analysis Point (ft)": sqrt((analysis_x - xw) ** 2 + (analysis_y - yw) ** 2),
        })

wheel_contribution_df = pd.DataFrame(wheel_contribution_rows)


col1, col2, col3, col4 = st.columns(4)

col1.metric("Aircraft", aircraft)
col2.metric("Main Gear Wheels", number_of_wheels)
col3.metric("Individual Tire Load", f"{wheel_load:,.0f} lb")
col4.metric("Calculated Contact Pressure", f"{contact_pressure_psi:,.1f} psi")


contact_radius_in = contact_radius_ft * 12.0
contact_diameter_in = contact_radius_in * 2.0

st.subheader("Model Assumptions")

st.write(
    f"Each aircraft tire is modeled as a **finite circular loaded area**.  \n"
    f"Contact radius = **{contact_radius_in:.1f} in**  \n"
    f"Contact diameter = **{contact_diameter_in:.1f} in**  \n"
    f"Contact area = **{contact_area_sq_in:,.1f} in²**  \n"
    f"Individual tire load = **{wheel_load:,.0f} lb**  \n"
    f"Calculated contact pressure = **{contact_pressure_psi:,.1f} psi**  \n"
    f"Analysis location = **({analysis_x:.2f}, {analysis_y:.2f}) ft**"
)

st.info(
    "In this version, the tire contact radius is fixed and contact pressure is calculated as "
    "individual wheel load divided by contact area. Therefore, increasing aircraft load increases "
    "calculated contact pressure and increases the stress response."
)


st.subheader("Single Wheel vs Combined Gear Stress")

fig_compare = px.line(
    comparison_df,
    x="Depth Below Pavement (ft)",
    y="Vertical Stress (psf)",
    color="Stress Type",
    markers=True,
    title=f"Stress Decay with Depth - {aircraft}",
)

fig_compare.update_layout(
    xaxis_title="Depth Below Pavement (ft)",
    yaxis_title="Vertical Stress (psf)",
)

st.plotly_chart(fig_compare, use_container_width=True)


st.subheader("Individual Wheel Contributions at Selected Analysis Point")

fig_wheels = px.line(
    wheel_contribution_df,
    x="Depth Below Pavement (ft)",
    y="Stress Contribution (psf)",
    color="Wheel",
    title="Individual Wheel Stress Contributions",
)

fig_wheels.update_layout(
    xaxis_title="Depth Below Pavement (ft)",
    yaxis_title="Stress Contribution (psf)",
)

st.plotly_chart(fig_wheels, use_container_width=True)


st.subheader("Plan-View Stress Bulb Overlap Map")

selected_depth = st.slider(
    "Selected Depth for Plan-View Stress Map (ft)",
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

x_grid = np.linspace(-grid_range, grid_range, 90)
y_grid = np.linspace(-grid_range, grid_range, 90)

grid_rows = []

for x in x_grid:
    for y in y_grid:
        total = total_stress_at_point(
            wheel_coordinates,
            contact_pressure_psf,
            contact_radius_ft,
            x,
            y,
            selected_depth,
        )

        grid_rows.append({
            "X (ft)": x,
            "Y (ft)": y,
            "Combined Stress (psf)": total,
        })

grid_df = pd.DataFrame(grid_rows)

fig_heat = px.density_heatmap(
    grid_df,
    x="X (ft)",
    y="Y (ft)",
    z="Combined Stress (psf)",
    nbinsx=90,
    nbinsy=90,
    title=f"Combined Stress Bulb Overlap at {selected_depth:.1f} ft Depth",
)

wheel_df = pd.DataFrame(wheel_coordinates, columns=["X (ft)", "Y (ft)"])

fig_wheel_points = px.scatter(
    wheel_df,
    x="X (ft)",
    y="Y (ft)",
)

for trace in fig_wheel_points.data:
    trace.name = "Wheel Locations"
    trace.marker.size = 10
    fig_heat.add_trace(trace)

fig_heat.add_scatter(
    x=[analysis_x],
    y=[analysis_y],
    mode="markers",
    marker=dict(size=14, symbol="x"),
    name="Analysis Point",
)

fig_heat.update_yaxes(scaleanchor="x", scaleratio=1)

st.plotly_chart(fig_heat, use_container_width=True)


st.subheader("Simplified Main Gear Layout")

fig_layout = px.scatter(
    wheel_df,
    x="X (ft)",
    y="Y (ft)",
    title=f"Simplified Main Gear Wheel Layout - {aircraft}",
)

fig_layout.add_scatter(
    x=[analysis_x],
    y=[analysis_y],
    mode="markers",
    marker=dict(size=14, symbol="x"),
    name="Analysis Point",
)

fig_layout.update_yaxes(scaleanchor="x", scaleratio=1)

st.plotly_chart(fig_layout, use_container_width=True)


st.subheader("Stress Table")

table_df = comparison_df.pivot(
    index="Depth Below Pavement (ft)",
    columns="Stress Type",
    values="Vertical Stress (psf)"
).reset_index()

st.dataframe(
    table_df.style.format({
        "Depth Below Pavement (ft)": "{:.2f}",
        "Single Wheel Centerline": "{:,.1f}",
        "Combined Gear at Selected Point": "{:,.1f}",
    }),
    use_container_width=True,
)


st.subheader("Wheel Contribution Table")

st.dataframe(
    wheel_contribution_df.style.format({
        "Depth Below Pavement (ft)": "{:.2f}",
        "Stress Contribution (psf)": "{:,.1f}",
        "Stress Contribution (psi)": "{:,.2f}",
        "Distance From Analysis Point (ft)": "{:.2f}",
    }),
    use_container_width=True,
)


st.warning(
    "Engineering note: This is an educational approximation. It uses fixed equivalent circular tire contact areas "
    "and superposes individual wheel stress bulbs. It does not replace layered elastic pavement analysis, "
    "FAARFIELD, finite element modeling, aircraft-specific tire-pavement contact mechanics, gear wander, "
    "load repetitions, or pass-to-coverage analysis."
)

