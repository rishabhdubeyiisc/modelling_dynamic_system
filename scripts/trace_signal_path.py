#!/usr/bin/env python3
"""
Complete State Variable Tracer for Generator Response
Traces ALL state variables to see the complete system response
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
    df = pd.read_csv(fname)  # Now reads proper CSV format instead of space-separated
    t = df["time"].values
    
    # Voltage control variables (handle missing columns)
    if "Vref" in df.columns:
        vref = df["Vref"].values
    else:
        vref = np.ones_like(t) * 1.04  # Estimate as 1.04 pu
    
    efd = df["Efd"].values
    eq1 = df["Eq1"].values      # Eq' (transient EMF)
    eq2 = df["Eq2"].values      # Eq'' (sub-transient EMF)
    ed1 = df["Ed1"].values      # Ed' (transient EMF)
    ed2 = df["Ed2"].values      # Ed'' (sub-transient EMF)
    vt = df["Vt"].values
    delta = df["delta"].values  # Rotor angle
    slip = df["slip"].values    # Rotor slip
    iq = df["iq_0"].values      # Q-axis current
    id = df["id_0"].values      # D-axis current
    
    # AVR error signal 
    if "EROR" in df.columns:
        error = df["EROR"].values
    else:
        error = vref - vt  # Calculate as Vref - Vt
    
    mech_power = df["mech_power"].values
    elec_power = df["elec_power"].values # Electrical power
    
    # Bus voltages (use the new column names)
    if "VQ_0" in df.columns:
        vq = df["VQ_0"].values
        vd = df["VD_0"].values
    else:
        vq = df["vq"].values      # Bus Q voltage
        vd = df["vd"].values      # Bus D voltage
    
    # Find the step location
    vref_diff = np.diff(vref)
    step_idx = np.argmax(np.abs(vref_diff)) + 1
    step_time = t[step_idx]
    
    print("==== COMPLETE STATE VARIABLE ANALYSIS ====")
    print(f"Step detected at t = {step_time:.2f} s")
    print()
    
    # Pre/post analysis
    pre_idx = max(0, step_idx - 50)
    post_idx = min(len(t) - 1, step_idx + 200)
    
    # Calculate pre/post values for ALL variables
    variables = {
        'Vref': vref, 'Efd': efd, 'Eq1': eq1, 'Eq2': eq2, 
        'Ed1': ed1, 'Ed2': ed2, 'Vt': vt, 'delta': delta,
        'slip': slip, 'iq': iq, 'id': id, 'P_elec': elec_power,
        'VQ': vq, 'VD': vd
    }
    
    print("Variable    Pre-Step  →  Post-Step   Change     % Change   Status")
    print("=" * 75)
    
    for name, var in variables.items():
        pre_val = np.mean(var[pre_idx:step_idx])
        post_val = np.mean(var[step_idx+10:step_idx+50])
        change = post_val - pre_val
        pct_change = 100 * change / pre_val if abs(pre_val) > 1e-10 else 0
        
        # Status assessment
        if abs(pct_change) > 50:
            status = "🟢 LARGE"
        elif abs(pct_change) > 10:
            status = "🟡 MEDIUM"
        elif abs(pct_change) > 1:
            status = "🟠 SMALL"
        else:
            status = "🔴 MINIMAL"
            
        print(f"{name:10s} {pre_val:8.3f}  →  {post_val:8.3f}   {change:+7.3f}   {pct_change:+6.1f}%   {status}")
    
    print()
    
    # Create comprehensive plot with ALL variables
    fig, axes = plt.subplots(4, 4, figsize=(20, 12))
    axes = axes.flatten()
    
    # Time window for plotting
    time_window = (t >= step_time - 5) & (t <= step_time + 15)
    t_plot = t[time_window]
    
    # Plot each variable
    plot_configs = [
        ('Vref', vref, 'blue', 'Vref (pu)'),
        ('Efd', efd, 'green', 'Efd (pu)'),
        ('Eq1 (Eq\')', eq1, 'magenta', 'Eq\' (pu)'),
        ('Eq2 (Eq\'\')', eq2, 'cyan', 'Eq\'\' (pu)'),
        ('Ed1 (Ed\')', ed1, 'orange', 'Ed\' (pu)'),
        ('Ed2 (Ed\'\')', ed2, 'red', 'Ed\'\' (pu)'),
        ('Vt', vt, 'black', 'Vt (pu)'),
        ('delta', delta, 'purple', 'δ (rad)'),
        ('slip', slip, 'olive', 'slip (pu)'),
        ('iq', iq, 'pink', 'iq (pu)'),
        ('id', id, 'gray', 'id (pu)'),
        ('P_elec', elec_power, 'brown', 'P_elec (pu)'),
        ('VQ', vq, 'pink', 'VQ (pu)'),
        ('VD', vd, 'gray', 'VD (pu)')
    ]
    
    for i, (title, var, color, ylabel) in enumerate(plot_configs):
        ax = axes[i]
        ax.plot(t_plot, var[time_window], color=color, linewidth=2)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.axvline(step_time, color='red', linestyle='--', alpha=0.7)
        ax.set_title(title, fontsize=10)
        
        # Add change annotation
        pre_val = np.mean(var[pre_idx:step_idx])
        post_val = np.mean(var[step_idx+10:step_idx+50])
        change = post_val - pre_val
        pct_change = 100 * change / pre_val if abs(pre_val) > 1e-10 else 0
        ax.text(0.02, 0.98, f'Δ{pct_change:+.1f}%', transform=ax.transAxes, 
                verticalalignment='top', fontsize=8, 
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Set x-label for bottom row and turn off unused subplots
    for i in range(len(plot_configs), len(axes)):
        axes[i].set_visible(False)  # Hide unused subplots
    
    for i in range(12, 16):
        if i < len(plot_configs):
            axes[i].set_xlabel('Time (s)')
    
    plt.suptitle('Complete State Variable Response to Vref Step', fontsize=16)
    plt.tight_layout()
    
    # Save plot
    plot_path = "report_plot/complete_state_trace.png"
    os.makedirs("report_plot", exist_ok=True)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    print(f"Complete state trace plot saved to: {plot_path}")
    
    # Signal chain health check
    print("\n==== SIGNAL CHAIN HEALTH CHECK ====")
    efd_change = abs(100 * (np.mean(efd[step_idx+10:step_idx+50]) - np.mean(efd[pre_idx:step_idx])) / np.mean(efd[pre_idx:step_idx]))
    eq1_change = abs(100 * (np.mean(eq1[step_idx+10:step_idx+50]) - np.mean(eq1[pre_idx:step_idx])) / np.mean(eq1[pre_idx:step_idx]))
    eq2_change = abs(100 * (np.mean(eq2[step_idx+10:step_idx+50]) - np.mean(eq2[pre_idx:step_idx])) / np.mean(eq2[pre_idx:step_idx]))
    vt_change = abs(100 * (np.mean(vt[step_idx+10:step_idx+50]) - np.mean(vt[pre_idx:step_idx])) / np.mean(vt[pre_idx:step_idx]))
    
    print(f"AVR → Efd:     {efd_change:6.1f}% {'✅ WORKING' if efd_change > 50 else '❌ BROKEN'}")
    print(f"Efd → Eq':     {eq1_change:6.1f}% {'✅ WORKING' if eq1_change > 10 else '❌ BROKEN'}")
    print(f"Eq' → Eq'':    {eq2_change:6.1f}% {'✅ WORKING' if eq2_change > 5 else '❌ BROKEN'}")
    print(f"Eq'' → Vt:     {vt_change:6.1f}% {'✅ WORKING' if vt_change > 5 else '❌ BROKEN'}")
    
    if efd_change > 50 and eq1_change < 10:
        print("\n🔍 DIAGNOSIS: Problem is in Efd → Eq' transfer (differential equation)")
    elif eq1_change > 10 and eq2_change < 5:
        print("\n🔍 DIAGNOSIS: Problem is in Eq' → Eq'' transfer (differential equation)")
    elif eq2_change > 5 and vt_change < 5:
        print("\n🔍 DIAGNOSIS: Problem is in Eq'' → Vt transfer (Norton current calculation)")
    elif efd_change < 50:
        print("\n🔍 DIAGNOSIS: Problem is in AVR (Vref → Efd)")
    else:
        print("\n✅ All signal chain stages appear to be working")

if __name__ == "__main__":
    main() 