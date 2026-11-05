import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
from crude_oil_fouling_interactive import simulate, RATE_SCALE_DEFAULT, AGEING_RATES

st.set_page_config(page_title="Crude-oil fouling simulator", layout="wide")
st.title("Crude-oil fouling: interactive simulation")
st.caption(
    "Reduced-order research prototype based on Coletti & Macchietto (2011). "
    "Scenario-based recalculation only; not a live plant-data connection."
)

data_path = Path(__file__).with_name("Figure8_digitised_daily_inputs.csv")

with st.sidebar:
    st.header("Simulation controls")
    days = st.slider("Simulation duration (days)", 30, 347, 347, 1)
    nodes = st.select_slider("Axial nodes per pass", options=[7, 13, 25, 49], value=13)
    ageing = st.selectbox("Deposit ageing scenario", list(AGEING_RATES), index=1)
    rate = st.slider("Fouling-rate multiplier", 0.01, 0.50, float(RATE_SCALE_DEFAULT), 0.005)
    shell_h = st.slider("Shell-side heat-transfer coefficient (W/m²K)", 100, 1000, 250, 25)
    tube_flow = st.slider("Tube flow multiplier", 0.5, 1.5, 1.0, 0.05)
    shell_flow = st.slider("Shell flow multiplier", 0.5, 1.5, 1.0, 0.05)
    run = st.button("Run / update simulation", type="primary")

if "results" not in st.session_state or run:
    with st.spinner("Solving the model..."):
        summary, local = simulate(
            data_path, days, nodes, ageing, rate, shell_h, tube_flow, shell_flow
        )
        st.session_state.results = (summary, local)

summary, local = st.session_state.results
last = summary.iloc[-1]

a, b, c, d = st.columns(4)
a.metric("Mean deposit thickness", f"{last.mean_deposit_thickness_mm:.3f} mm")
b.metric("Maximum deposit thickness", f"{last.max_deposit_thickness_mm:.3f} mm")
c.metric("Raw pressure-drop estimate", f"{last.raw_pressure_drop_bar:.3f} bar")
d.metric("Modelled tube outlet", f"{last.tube_outlet_model_C:.1f} °C")

st.warning(
    "Outputs are provisional model estimates. Input profiles are digitised from a "
    "published figure; several property and heat-transfer assumptions remain "
    "approximate. Energy closure and physical validation must be reviewed before "
    "using results as engineering predictions."
)

tab1, tab2, tab3 = st.tabs(["Time trends", "Axial profiles", "Download data"])
with tab1:
    st.plotly_chart(
        px.line(
            summary,
            x="day",
            y=["mean_deposit_thickness_mm", "max_deposit_thickness_mm"],
            labels={"value": "Thickness (mm)", "variable": "Metric"},
            title="Deposit thickness vs time",
        ),
        use_container_width=True,
    )
    st.plotly_chart(
        px.line(
            summary,
            x="day",
            y="raw_pressure_drop_bar",
            title="Raw estimated pressure drop vs time",
            labels={"raw_pressure_drop_bar": "Pressure drop (bar)"},
        ),
        use_container_width=True,
    )
    st.plotly_chart(
        px.line(
            summary,
            x="day",
            y=["tube_duty_MW", "shell_duty_from_inputs_MW"],
            title="Heat-duty comparison (diagnostic)",
            labels={"value": "Duty (MW)", "variable": "Side"},
        ),
        use_container_width=True,
    )

with tab2:
    selected_day = int(summary.day.max())
    profile = local[local.day == selected_day]
    st.plotly_chart(
        px.line(
            profile,
            x="z_m",
            y="deposit_m",
            color="pass",
            facet_col="shell",
            markers=True,
            title=f"Axial deposit profile at day {selected_day}",
            labels={"deposit_m": "Deposit thickness (m)", "z_m": "Axial coordinate (m)"},
        ),
        use_container_width=True,
    )
    st.plotly_chart(
        px.line(
            profile,
            x="z_m",
            y="film_C",
            color="pass",
            facet_col="shell",
            title="Local film-temperature profile",
        ),
        use_container_width=True,
    )

with tab3:
    st.download_button(
        "Download daily summary CSV",
        summary.to_csv(index=False).encode("utf-8"),
        "daily_summary.csv",
        "text/csv",
    )
    st.download_button(
        "Download local axial results CSV",
        local.to_csv(index=False).encode("utf-8"),
        "local_axial_results.csv",
        "text/csv",
    )
    st.dataframe(summary.tail(10), use_container_width=True)
