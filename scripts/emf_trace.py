#!/usr/bin/env python3
"""
EMF Variable Tracer
Tracks Eq1, Eq2, Ed1, Ed2 and their differential equations to find where the update breaks
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

plt.switch_backend("Agg")

def main():
    fname = "PLOT_SAVED_BUS.csv"
    if not os.path.exists(fname):
        raise SystemExit(f"❌ {fname} not found – run the simulation first.")

    # Read data
    df = pd.read_csv(fname)
    t = df["time"].values
    
    # EMF Variables
    efd = df["Efd"].values
    eq1 = df["Eq1"].values      # Eq' (transient EMF)
    eq2 = df["Eq2"].values      # Eq'' (sub-transient EMF)
    ed1 = df["Ed1"].values      # Ed' (transient EMF)
    ed2 = df["Ed2"].values      # Ed'' (sub-transient EMF)
    
    # Currents
    iq = df["iq_0"].values
    id = df["id_0"].values
    
    # Voltage control variables (estimate Vref if not available)
    if "Vref" in df.columns:
        vref = df["Vref"].values
    else:
        # Estimate Vref as constant 1.04 pu (typical generator reference)
        vref = np.ones_like(t) * 1.04
    
    vt = df["Vt"].values
    
    if "EROR" in df.columns:
        error = df["EROR"].values
    else:
        # Calculate voltage error as Vref - Vt
        error = vref - vt
    
    # Find the step location
    vref_diff = np.diff(vref)
    step_idx = np.argmax(np.abs(vref_diff)) + 1
    step_time = t[step_idx]
    
    print("==== EMF VARIABLE DETAILED TRACE ====")
    print(f"Step detected at t = {step_time:.2f} s")
    print()
    
    # Show values before and after step
    pre_idx = max(0, step_idx - 10)
    post_idx = min(len(t) - 1, step_idx + 50)
    
    print("TIME SERIES AROUND STEP:")
    print("Time    Vref     Efd      Eq1      Eq2      Ed1      Ed2      iq       id")
    print("=" * 85)
    
    for i in range(max(0, step_idx-5), min(len(t), step_idx+15)):
        marker = " >>> " if i == step_idx else "     "
        print(f"{t[i]:6.1f}{marker}{vref[i]:7.3f}  {efd[i]:7.3f}  {eq1[i]:7.3f}  {eq2[i]:7.3f}  {ed1[i]:7.3f}  {ed2[i]:7.3f}  {iq[i]:7.3f}  {id[i]:7.3f}")
    
    print()
    print("==== DIFFERENTIAL EQUATION ANALYSIS ====")
    
    # Calculate theoretical derivatives
    # From the differential equations in F_VECTOR_CALC_single:
    # dEq'/dt = (1/Tdo1) * (Efd - Eq' + id*(Xd0-Xd1))
    # dEd'/dt = (-1/Tqo1) * (Ed' + iq*(Xq0-Xq1))
    # dEq''/dt = (1/Tdo2) * (Eq' - Eq'' + id*(Xd1-Xd2))
    # dEd''/dt = (1/Tqo2) * (Ed' - Ed'' + iq*(Xq2-Xq1))
    
    # Generator parameters (from System_Data/system_data.txt line 2)
    # 23.64 0.146 0.100 0.060 0.0969 0.085 0.070 1.0 1.0 0.31 0.31
    Xd0 = 0.146
    Xd1 = 0.100  
    Xd2 = 0.060
    Xq0 = 0.0969
    Xq1 = 0.085
    Xq2 = 0.070
    Tdo1 = 1.0
    Tqo1 = 1.0
    Tdo2 = 0.31
    Tqo2 = 0.31
    
    # Calculate derivatives at step point
    step_efd = efd[step_idx]
    step_eq1 = eq1[step_idx]
    step_eq2 = eq2[step_idx]
    step_ed1 = ed1[step_idx]
    step_ed2 = ed2[step_idx]
    step_id = id[step_idx]
    step_iq = iq[step_idx]
    
    # Theoretical derivatives
    dEq1_dt_theory = (1/Tdo1) * (step_efd - step_eq1 + step_id*(Xd0-Xd1))
    dEd1_dt_theory = (-1/Tqo1) * (step_ed1 + step_iq*(Xq0-Xq1))
    dEq2_dt_theory = (1/Tdo2) * (step_eq1 - step_eq2 + step_id*(Xd1-Xd2))
    dEd2_dt_theory = (1/Tqo2) * (step_ed1 - step_ed2 + step_iq*(Xq2-Xq1))
    
    print(f"At step time t = {step_time:.2f} s:")
    print(f"Efd = {step_efd:.3f}")
    print(f"Eq1 = {step_eq1:.3f}")
    print(f"Eq2 = {step_eq2:.3f}")
    print(f"Ed1 = {step_ed1:.3f}")
    print(f"Ed2 = {step_ed2:.3f}")
    print(f"id  = {step_id:.3f}")
    print(f"iq  = {step_iq:.3f}")
    print()
    
    print("THEORETICAL DERIVATIVES:")
    print(f"dEq'/dt  = (1/{Tdo1}) * ({step_efd:.3f} - {step_eq1:.3f} + {step_id:.3f}*{Xd0-Xd1:.3f}) = {dEq1_dt_theory:.6f}")
    print(f"dEd'/dt  = (-1/{Tqo1}) * ({step_ed1:.3f} + {step_iq:.3f}*{Xq0-Xq1:.3f}) = {dEd1_dt_theory:.6f}")
    print(f"dEq''/dt = (1/{Tdo2}) * ({step_eq1:.3f} - {step_eq2:.3f} + {step_id:.3f}*{Xd1-Xd2:.3f}) = {dEq2_dt_theory:.6f}")
    print(f"dEd''/dt = (1/{Tqo2}) * ({step_ed1:.3f} - {step_ed2:.3f} + {step_iq:.3f}*{Xq2-Xq1:.3f}) = {dEd2_dt_theory:.6f}")
    print()
    
    # Calculate actual derivatives from data
    dt = t[1] - t[0]  # time step
    if step_idx < len(t) - 1:
        dEq1_dt_actual = (eq1[step_idx+1] - eq1[step_idx]) / dt
        dEd1_dt_actual = (ed1[step_idx+1] - ed1[step_idx]) / dt
        dEq2_dt_actual = (eq2[step_idx+1] - eq2[step_idx]) / dt
        dEd2_dt_actual = (ed2[step_idx+1] - ed2[step_idx]) / dt
        
        print("ACTUAL DERIVATIVES FROM DATA:")
        print(f"dEq'/dt  = {dEq1_dt_actual:.6f}")
        print(f"dEd'/dt  = {dEd1_dt_actual:.6f}")
        print(f"dEq''/dt = {dEq2_dt_actual:.6f}")
        print(f"dEd''/dt = {dEd2_dt_actual:.6f}")
        print()
        
        print("DERIVATIVE COMPARISON:")
        print(f"dEq'/dt:  Theory = {dEq1_dt_theory:.6f}, Actual = {dEq1_dt_actual:.6f}, Ratio = {dEq1_dt_actual/dEq1_dt_theory if abs(dEq1_dt_theory) > 1e-10 else 'N/A'}")
        print(f"dEd'/dt:  Theory = {dEd1_dt_theory:.6f}, Actual = {dEd1_dt_actual:.6f}, Ratio = {dEd1_dt_actual/dEd1_dt_theory if abs(dEd1_dt_theory) > 1e-10 else 'N/A'}")
        print(f"dEq''/dt: Theory = {dEq2_dt_theory:.6f}, Actual = {dEq2_dt_actual:.6f}, Ratio = {dEq2_dt_actual/dEq2_dt_theory if abs(dEq2_dt_theory) > 1e-10 else 'N/A'}")
        print(f"dEd''/dt: Theory = {dEd2_dt_theory:.6f}, Actual = {dEd2_dt_actual:.6f}, Ratio = {dEd2_dt_actual/dEd2_dt_theory if abs(dEd2_dt_theory) > 1e-10 else 'N/A'}")
        print()
    
    # Check if derivatives are effectively zero
    derivative_threshold = 1e-6
    print("DERIVATIVE MAGNITUDE CHECK:")
    print(f"dEq'/dt magnitude: {abs(dEq1_dt_theory):.2e} {'❌ TOO SMALL' if abs(dEq1_dt_theory) < derivative_threshold else '✅ SIGNIFICANT'}")
    print(f"dEd'/dt magnitude: {abs(dEd1_dt_theory):.2e} {'❌ TOO SMALL' if abs(dEd1_dt_theory) < derivative_threshold else '✅ SIGNIFICANT'}")
    print(f"dEq''/dt magnitude: {abs(dEq2_dt_theory):.2e} {'❌ TOO SMALL' if abs(dEq2_dt_theory) < derivative_threshold else '✅ SIGNIFICANT'}")
    print(f"dEd''/dt magnitude: {abs(dEd2_dt_theory):.2e} {'❌ TOO SMALL' if abs(dEd2_dt_theory) < derivative_threshold else '✅ SIGNIFICANT'}")
    print()
    
    # Analyze the terms
    print("DETAILED TERM ANALYSIS FOR dEq'/dt:")
    print(f"  Main drive term: (Efd - Eq') = ({step_efd:.3f} - {step_eq1:.3f}) = {step_efd - step_eq1:.3f}")
    print(f"  Current term: id*(Xd0-Xd1) = {step_id:.3f}*{Xd0-Xd1:.3f} = {step_id*(Xd0-Xd1):.3f}")
    print(f"  Total numerator: {(step_efd - step_eq1 + step_id*(Xd0-Xd1)):.3f}")
    print(f"  Time constant factor: 1/Tdo1 = 1/{Tdo1} = {1/Tdo1:.3f}")
    print(f"  Final derivative: {dEq1_dt_theory:.6f}")
    
    if abs(dEq1_dt_theory) < derivative_threshold:
        print("  🚨 PROBLEM: Derivative is essentially zero!")
        if abs(step_efd - step_eq1) < 0.1:
            print("     - Efd and Eq1 are too close")
        if abs(step_id*(Xd0-Xd1)) > abs(step_efd - step_eq1):
            print("     - Current term dominates and cancels main term")
    
    print()
    
    # Create detailed plot
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    
    time_window = (t >= step_time - 5) & (t <= step_time + 20)
    t_plot = t[time_window]
    
    # Plot EMF variables
    axes[0,0].plot(t_plot, efd[time_window], 'red', linewidth=2, label='Efd')
    axes[0,0].set_ylabel('Efd (pu)')
    axes[0,0].set_title('Field EMF (Efd)')
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[0,0].legend()
    
    axes[0,1].plot(t_plot, eq1[time_window], 'blue', linewidth=2, label="Eq'")
    axes[0,1].plot(t_plot, eq2[time_window], 'cyan', linewidth=2, label='Eq"')
    axes[0,1].set_ylabel('EMF (pu)')
    axes[0,1].set_title('Q-axis EMFs')
    axes[0,1].grid(True, alpha=0.3)
    axes[0,1].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[0,1].legend()
    
    axes[1,0].plot(t_plot, ed1[time_window], 'orange', linewidth=2, label="Ed'")
    axes[1,0].plot(t_plot, ed2[time_window], 'red', linewidth=2, label='Ed"')
    axes[1,0].set_ylabel('EMF (pu)')
    axes[1,0].set_title('D-axis EMFs')
    axes[1,0].grid(True, alpha=0.3)
    axes[1,0].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[1,0].legend()
    
    axes[1,1].plot(t_plot, iq[time_window], 'purple', linewidth=2, label='iq')
    axes[1,1].plot(t_plot, id[time_window], 'brown', linewidth=2, label='id')
    axes[1,1].set_ylabel('Current (pu)')
    axes[1,1].set_title('Currents')
    axes[1,1].grid(True, alpha=0.3)
    axes[1,1].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[1,1].legend()
    
    axes[2,0].plot(t_plot, vref[time_window], 'green', linewidth=2, label='Vref')
    axes[2,0].plot(t_plot, vt[time_window], 'blue', linewidth=2, label='Vt')
    axes[2,0].set_ylabel('Voltage (pu)')
    axes[2,0].set_title('Reference vs Terminal Voltage')
    axes[2,0].set_xlabel('Time (s)')
    axes[2,0].grid(True, alpha=0.3)
    axes[2,0].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[2,0].legend()
    
    axes[2,1].plot(t_plot, error[time_window], 'red', linewidth=2, label='Error')
    axes[2,1].set_ylabel('Error (pu)')
    axes[2,1].set_title('AVR Error Signal')
    axes[2,1].set_xlabel('Time (s)')
    axes[2,1].grid(True, alpha=0.3)
    axes[2,1].axvline(step_time, color='black', linestyle='--', alpha=0.7)
    axes[2,1].legend()
    
    plt.suptitle('Complete EMF Variable Trace', fontsize=16)
    plt.tight_layout()
    
    # Save plot
    plot_path = "report_plot/emf_trace.png"
    os.makedirs("report_plot", exist_ok=True)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"EMF trace plot saved to: {plot_path}")

if __name__ == "__main__":
    main() 