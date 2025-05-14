import os
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pipeline_model import solve_pipeline

# ---------------------
# Page configuration
# ---------------------
st.set_page_config(
    page_title="Mixed Integer Nonlinear Convex Optimization of Pipeline Operations",
    layout="wide"
)

# ---------------------
# Custom CSS
# ---------------------
st.markdown(
    '''
    <style>
      .metric-card {
        background: #ffffff;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        height: 120px;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        text-align: center;
        white-space: nowrap;
      }
      .metric-card h3 { margin: 0; font-size: 16px; font-weight: 600; }
      .metric-card h1 { margin: 0; font-size: 24px; font-weight: bold; }
      .stTabs [role=tab] { font-weight: bold; }
    </style>
    ''', unsafe_allow_html=True
)

# ---------------------
# Header
# ---------------------
st.markdown(
    """
    <div style='background: linear-gradient(90deg, #4f81bd, #2c3e50); padding: 20px; border-radius:10px; margin-bottom: 20px;'>
      <h1 style='color: white; margin: 0; text-align: center;'>
        Mixed Integer Nonlinear Convex Optimization of Pipeline Operations
      </h1>
    </div>
    """, unsafe_allow_html=True
)

# ---------------------
# Sidebar Inputs
# ---------------------
with st.sidebar:
    st.header("Pipeline Inputs")
    with st.expander("Adjust Parameters", expanded=True):
        FLOW      = st.number_input("Flow rate (m3/hr)",      value=1000.0, min_value=0.0, step=0.1, format="%.2f")
        KV        = st.number_input("Kinematic Viscosity (cSt)", value=1.0,    min_value=0.0, step=0.01, format="%.2f")
        rho       = st.number_input("Density (kg/m3)",        value=850.0,  min_value=0.0, step=0.1, format="%.2f")
        SFC_J     = st.number_input("SFC Jamnagar (gm/bhp/hr)", value=210.0,  min_value=0.0, step=0.1, format="%.2f")
        SFC_R     = st.number_input("SFC Rajkot (gm/bhp/hr)",   value=215.0,  min_value=0.0, step=0.1, format="%.2f")
        SFC_S     = st.number_input("SFC Surendranagar (gm/bhp/hr)", value=220.0, min_value=0.0, step=0.1, format="%.2f")
        RateDRA   = st.number_input("DRA Rate (INR/L)",        value=1.0,    min_value=0.0, step=0.01, format="%.2f")
        Price_HSD = st.number_input("HSD Rate (INR/L)",        value=90.0,   min_value=0.0, step=0.1, format="%.2f")
	neos_email = st.text_input("NEOS Email", value="parichay.nitwarangal@gmail.com")
    	run = st.button("Run Optimization")



# before calling solver, ensure NEOS_EMAIL is set
if run:
    if not neos_email:
        st.error("Please enter your NEOS-registered email address.")
    else:
        os.environ['NEOS_EMAIL'] = neos_email

        with st.spinner("Submitting to NEOS… this can take up to 2–3 minutes"):
            try:
                res = solve_pipeline(FLOW, KV, rho, SFC_J, SFC_R, SFC_S, RateDRA, Price_HSD)
            except Exception as e:
                st.error(f"Solver failed: {e}")
                st.stop()
    os.environ['NEOS_EMAIL'] = neos_email

    with st.spinner("Solving optimization..."):
        res = solve_pipeline(…)


    # KPIs
    total = res.get('total_cost', 0)
    pumps_total = sum(int(res.get(f"num_pumps_{s.lower()}",0)) for s in stations)
    speeds = [res.get(f"speed_{s.lower()}",0) for s in stations if res.get(f"num_pumps_{s.lower()}",0)>0]
    effs   = [res.get(f"efficiency_{s.lower()}",0) for s in stations if res.get(f"num_pumps_{s.lower()}",0)>0]
    avg_speed = np.mean(speeds) if speeds else 0
    avg_eff   = np.mean(effs)   if effs else 0

    cols = st.columns(5)
    cols[0].markdown(f"<div class='metric-card'><h3>Total Cost</h3><h1>₹{total:,.2f}</h1></div>", unsafe_allow_html=True)
    cols[1].markdown(f"<div class='metric-card'><h3>Total Pumps</h3><h1>{pumps_total}</h1></div>", unsafe_allow_html=True)
    cols[2].markdown(f"<div class='metric-card'><h3>Avg Pump Speed</h3><h1>{avg_speed:.2f} rpm</h1></div>", unsafe_allow_html=True)
    cols[3].markdown(f"<div class='metric-card'><h3>Avg Pump Efficiency</h3><h1>{avg_eff:.2f}%</h1></div>", unsafe_allow_html=True)
    cols[4].markdown(f"<div class='metric-card'><h3>Flow Rate</h3><h1>{FLOW:.2f} m3/hr</h1></div>", unsafe_allow_html=True)

    # Tabs
    tab1, tab2, tab3 = st.tabs(["Summary Table","Cost Charts","Performance Charts"])

    # Build DataFrame
    data = {"Process particulars":[
        "Power & Fuel cost (INR/day)","DRA cost (INR/day)","Total Operating Cost (INR/day)",
        "No. of Pumps","Pump Speed (rpm)","Pump Efficiency (%)","Reynold's no.",
        "Dynamic Head Loss (mcl)","Velocity (m/s)","RH (mcl)","SDH (mcl)","Drag reduction (%)"
    ]}
    for s in stations:
        k = s.lower()
        pumps = int(res.get(f"num_pumps_{k}",0))
        speed = res.get(f"speed_{k}",0) if pumps>0 else 0
        eff = res.get(f"efficiency_{k}",0) if pumps>0 else 0
        data[s] = [
            res.get(f"power_cost_{k}",0) + res.get(f"dra_cost_{k}",0),
            res.get(f"dra_cost_{k}",0),
            total if s==stations[0] else '-',
            pumps, speed, eff,
            res.get(f"reynolds_{k}",0), res.get(f"head_loss_{k}",0), res.get(f"velocity_{k}",0),
            res.get(f"residual_head_{k}",0), res.get(f"sdh_{k}",0), res.get(f"drag_reduction_{k}",0)
        ]
    df = pd.DataFrame(data)

    # Formatters
    formatters = {}
    for col in df.columns:
        if col == "Process particulars": continue
        if col == "No. of Pumps":
            formatters[col] = lambda x: f"{int(x)}" if isinstance(x, (int, float, np.integer, np.floating)) else '-'
        else:
            formatters[col] = lambda x: f"{x:,.2f}" if isinstance(x, (int, float, np.integer, np.floating)) else '-'

    # Tab1: Summary Table
    with tab1:
        st.subheader("Optimized Parameters Summary")
        styled_df = df.style.format(formatters).set_properties(**{'text-align':'center'})
        st.dataframe(styled_df, use_container_width=True)
        st.download_button("Download CSV", df.to_csv(index=False).encode(), "pipeline_report.csv", "text/csv")

    # Tab2: Cost Charts
    with tab2:
        st.subheader("Cost Breakdown per Station")
        cost_df = pd.DataFrame({
            'Station':     stations,
            'Power & Fuel':[res.get(f"power_cost_{s.lower()}",0) for s in stations],
            'DRA':         [res.get(f"dra_cost_{s.lower()}",0) for s in stations]
        })
        fig = px.bar(
            cost_df.melt(id_vars='Station',value_vars=['Power & Fuel','DRA']),
            x='Station', y='value', color='variable', barmode='stack', height=500
        )
        st.plotly_chart(fig, use_container_width=True)

    # Tab3: Performance Charts
    with tab3:
        st.subheader("Performance Metrics Over Sections")
        perf_df = pd.DataFrame({
            'Station':    stations,
            'Head Loss':  [res.get(f"head_loss_{s.lower()}",0) for s in stations],
            'Speed':      [res.get(f"speed_{s.lower()}",0) for s in stations],
            'Efficiency': [res.get(f"efficiency_{s.lower()}",0) for s in stations]
        })
        fig2 = go.Figure()
        fig2.add_trace(go.Bar(x=perf_df['Station'],y=perf_df['Head Loss'],name='Head Loss (mcl)'))
        fig2.add_trace(go.Scatter(x=perf_df['Station'],y=perf_df['Speed'],name='Pump Speed (rpm)',yaxis='y2',mode='lines+markers'))
        fig2.add_trace(go.Scatter(x=perf_df['Station'],y=perf_df['Efficiency'],name='Pump Efficiency (%)',yaxis='y2',mode='lines+markers'))
        fig2.update_layout(
            yaxis2=dict(overlaying='y', side='right', title='Speed/Efficiency'),
            xaxis_title='Station',
            height=500
        )
        st.plotly_chart(fig2, use_container_width=True)

    st.markdown("---")
    st.caption("© 2025 Developed by Parichay Das. All rights reserved.")
