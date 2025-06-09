#!/usr/bin/env python3
"""
Vref-to-Voltage Analysis Tool
Reads PLOT_SAVED_BUS.csv and produces:
  • Time-series overlay of Vref, Vt and Efd with the Vref-step highlighted
  • Scatter plots Vref vs Vt and Vref vs Efd with linear-fit slope and R²
  • Console report of pre-/post-step averages, steady-state error and correlation coefficients.
Outputs figures in report_plot/.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

plt.switch_backend("Agg")

def main():
    fname = "PLOT_SAVED_BUS.csv"
    if not os.path.exists(fname):
        raise SystemExit("❌ {} not found – run the simulation first.".format(fname))

    df = pd.read_csv(fname)  # Now reads proper CSV format instead of space-separated
    t = df["time"].values
    if "Vref" in df.columns:
        vref = df["Vref"].values
    else:
        # Estimate Vref as constant 1.04 pu for generator
        vref = np.ones_like(t) * 1.04
    vt = df["Vt"].values
    efd = df["Efd"].values

    # Detect Vref step – first time derivative > threshold
    dvref = np.diff(vref)
    step_idx = np.argmax(np.abs(dvref) > 1e-3)  # threshold 0.001 pu
    step_time = t[step_idx+1]

    # Check if a real step was detected
    max_vref_change = np.max(np.abs(dvref))
    if max_vref_change < 1e-6:  # No significant step found
        print("📝 No Vref step detected - analyzing steady-state performance")
        # Use first and last quarters for comparison
        quarter_len = len(t) // 4
        pre_mask = slice(0, quarter_len)
        post_mask = slice(-quarter_len, None)
        step_time = t[len(t)//2]  # Use middle of simulation as nominal "step time"
        
        pre_vref, post_vref = np.mean(vref[pre_mask]), np.mean(vref[post_mask])
        pre_vt,   post_vt   = np.mean(vt[pre_mask]),   np.mean(vt[post_mask])
        pre_efd,  post_efd  = np.mean(efd[pre_mask]),  np.mean(efd[post_mask])
        
    else:
        # Real step detected - use time-based masks
        pre_mask = t < step_time - 0.5
        post_mask = t > step_time + 30.0   # allow settling 30 s
        
        # Ensure masks are not empty
        if np.sum(pre_mask) == 0:
            pre_mask = t < step_time
        if np.sum(post_mask) == 0:
            post_mask = t > step_time
            
        # Handle case where masks might still be empty
        if np.sum(pre_mask) > 0 and np.sum(post_mask) > 0:
            pre_vref, post_vref = np.mean(vref[pre_mask]), np.mean(vref[post_mask])
            pre_vt,   post_vt   = np.mean(vt[pre_mask]),   np.mean(vt[post_mask])
            pre_efd,  post_efd  = np.mean(efd[pre_mask]),  np.mean(efd[post_mask])
        else:
            # Fallback: use first/last values
            pre_vref, post_vref = vref[0], vref[-1]
            pre_vt,   post_vt   = vt[0],   vt[-1]
            pre_efd,  post_efd  = efd[0],  efd[-1]

    # Correlation & regression (handle constant Vref case)
    vref_std = np.std(vref)
    if vref_std < 1e-6:  # Vref is essentially constant
        print("📝 Note: Vref is constant - analyzing voltage regulation performance instead")
        slope_vt, r_vt = float('nan'), float('nan')
        slope_efd, r_efd = float('nan'), float('nan')
        
        # Calculate voltage regulation metrics instead
        vt_std = np.std(vt)
        efd_std = np.std(efd)
        vt_regulation_error = np.max(np.abs(vt - np.mean(vref)))
        
        print(f"Voltage regulation analysis:")
        print(f"  Vref (constant): {np.mean(vref):.4f} pu")
        print(f"  Vt std deviation: {vt_std:.6f} pu")
        print(f"  Efd std deviation: {efd_std:.6f} pu") 
        print(f"  Max Vt error from Vref: {vt_regulation_error:.6f} pu")
        
    else:
        # Normal case: Vref varies, can do regression
        slope_vt, _, r_vt, _, _  = linregress(vref, vt)
        slope_efd, _, r_efd, _, _ = linregress(vref, efd)

    print("==== Vref Analysis Report ====")
    if max_vref_change > 1e-6:
        print(f"Step detected at t = {step_time:.2f} s (index {step_idx})")
    else:
        print(f"No step detected - steady-state analysis from t = {step_time:.2f} s")
    print(f"Vref : {pre_vref:.4f} → {post_vref:.4f}  (Δ = {post_vref-pre_vref:.4f})")
    print(f"Vt   : {pre_vt:.4f} → {post_vt:.4f}    (Δ = {post_vt-pre_vt:.4f})")
    print(f"Efd  : {pre_efd:.4f} → {post_efd:.4f}   (Δ = {post_efd-pre_efd:.4f})")
    
    if not np.isnan(slope_vt):
        print(f"Linear fit Vt = {slope_vt:.3f}*Vref  |  R = {r_vt:.4f}")
        print(f"Linear fit Efd = {slope_efd:.3f}*Vref |  R = {r_efd:.4f}")
    else:
        print("Linear regression skipped (Vref is constant)")

    # Plot time-series
    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(t, vref, label="Vref", lw=1.5)
    ax.plot(t, vt, label="Vt", lw=1.5)
    ax.plot(t, efd, label="Efd", lw=1.2)
    ax.axvline(step_time, color="k", ls="--", alpha=0.7, label="Vref step")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Voltage (pu)")
    ax.set_title("Vref Step Response")
    ax.grid(alpha=0.3)
    ax.legend()
    os.makedirs("report_plot", exist_ok=True)
    fig.savefig("report_plot/vref_step_response.png", dpi=300, bbox_inches="tight")

    # Scatter plots
    fig2, axes = plt.subplots(1,2, figsize=(12,5))
    axes[0].scatter(vref, vt, s=4, alpha=0.6)
    axes[0].set_xlabel("Vref (pu)")
    axes[0].set_ylabel("Vt (pu)")
    if not np.isnan(slope_vt):
        axes[0].set_title(f"Vt vs Vref  (slope={slope_vt:.3f}, R={r_vt:.3f})")
    else:
        axes[0].set_title("Vt vs Vref  (Vref constant)")
    axes[0].grid(alpha=0.3)

    axes[1].scatter(vref, efd, s=4, alpha=0.6, color="orange")
    axes[1].set_xlabel("Vref (pu)")
    axes[1].set_ylabel("Efd (pu)")
    if not np.isnan(slope_efd):
        axes[1].set_title(f"Efd vs Vref  (slope={slope_efd:.3f}, R={r_efd:.3f})")
    else:
        axes[1].set_title("Efd vs Vref  (Vref constant)")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    fig2.savefig("report_plot/vref_scatter.png", dpi=300, bbox_inches="tight")
    print("Plots saved to report_plot/: vref_step_response.png, vref_scatter.png")

if __name__ == "__main__":
    main() 