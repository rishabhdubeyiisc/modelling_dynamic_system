#!/usr/bin/env python3
"""
Power System Dynamics - Comprehensive Plotting Suite
Combines all plotting functionality for simulation analysis
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import os

# Set matplotlib to use non-interactive backend
plt.switch_backend('Agg')

# ------------------------------------------------------------------
# Helper: robust CSV loader that works with both PLOT.csv and
# PLOT_SAVED_BUS.csv produced by the simulator. Handles space-separated
# files without headers and cleans non-numeric rows.
# ------------------------------------------------------------------

def load_plot_csv():
    candidate_names = [
        'PLOT_SAVED_BUS.csv',
        'PLOT.csv',
    ]

    for fname in candidate_names:
        if os.path.exists(fname):
            try:
                df = pd.read_csv(
                    fname,
                    sep=',',
                    header=0
                )

                # Convert to numeric, coercing errors to NaN
                for col in df.columns:
                    if col != 'time':  # Keep time as is
                        df[col] = pd.to_numeric(df[col], errors='coerce')

                # Drop rows with any NaN values (bad data)
                df.dropna(inplace=True)
                df.reset_index(drop=True, inplace=True)

                print(f"Loaded data from '{fname}' with shape {df.shape}")
                return df
            except Exception as e:
                print(f"Failed reading {fname}: {e}")
                continue

    raise FileNotFoundError("Neither PLOT_SAVED_BUS.csv nor PLOT.csv found")

def _read_fault_parameters_from_main():
    """Parse main.c to obtain fault parameters; return defaults if not found."""
    fault_cycles = 8.0  # Default
    clear_disturbances = 150.0  # Default in seconds
    
    try:
        with open('main.c', 'r', encoding='utf-8') as f:
            text = f.read()
        
        # Read FAULT_CYCLES
        m = re.search(r"FAULT_CYCLES\s*=\s*([0-9]*\.?[0-9]+)", text)
        if m:
            fault_cycles = float(m.group(1))
        
        # Read CLEAR_DISTURBANCES (should be in seconds)
        m = re.search(r"CLEAR_DISTURBANCES\s*=\s*([0-9]*\.?[0-9]+)\s*\*\s*60", text)
        if m:
            clear_disturbances = float(m.group(1)) * 60  # Convert minutes to seconds
        else:
            # Try direct seconds value
            m = re.search(r"CLEAR_DISTURBANCES\s*=\s*([0-9]*\.?[0-9]+)", text)
            if m:
                clear_disturbances = float(m.group(1))
                
    except (FileNotFoundError, IOError):
        pass
    
    return fault_cycles, clear_disturbances

# New: Event log helpers
from psd_utils import get_event_times

# ------------------------------------------------------------------
# Fault timing detection now prefers explicit entries from sim/events.csv.
# If the events log is missing or incomplete, we gracefully fall back to the
# previous heuristic/data-driven approach.
# ------------------------------------------------------------------

def _detect_fault_timing_from_data(df, time):
    """Return (fault_start, fault_end) times.

    Priority order:
    1. Use explicit times from sim/events.csv (fast & reliable).
    2. Fallback to heuristic voltage-drop detection around code-scheduled
       timings (legacy behaviour).
    """

    # 1) Try direct event log first
    evt_start = get_event_times('fault_start')
    evt_end   = get_event_times('fault_end')
    if evt_start and evt_end:
        return evt_start[0], evt_end[0]
    
    # 2) Fallback to heuristic voltage-drop detection around code-scheduled timings
    fault_cycles, clear_disturbances = _read_fault_parameters_from_main()
    code_fault_start = 2 * clear_disturbances
    code_fault_end = code_fault_start + (fault_cycles / 60.0)
    
    # Look for significant voltage drops around the expected time
    vt = df['Vt'].values
    
    # Only search in a window around the expected fault time (±50 seconds)
    search_start_time = max(0, code_fault_start - 50)
    search_end_time = min(time[-1], code_fault_end + 50)
    
    search_mask = (time >= search_start_time) & (time <= search_end_time)
    search_indices = np.where(search_mask)[0]
    
    if len(search_indices) > 1:
        search_vt = vt[search_indices]
        search_time = time[search_indices]
        
        # Find the largest voltage drop in this window
        voltage_diff = np.diff(search_vt)
        significant_drops = np.where(voltage_diff < -0.1)[0]  # Drops > 0.1 pu
        
        if len(significant_drops) > 0:
            # Use the drop closest to expected fault time
            expected_fault_idx = np.argmin(np.abs(search_time - code_fault_start))
            drop_distances = np.abs(significant_drops - expected_fault_idx)
            best_drop_idx = significant_drops[np.argmin(drop_distances)]
            
            fault_start_idx = search_indices[best_drop_idx]
            fault_start_time = time[fault_start_idx]
            
            # Find recovery point (voltage starts increasing again)
            recovery_search = voltage_diff[best_drop_idx:]
            recovery_points = np.where(recovery_search > 0.05)[0]
            if len(recovery_points) > 0:
                fault_end_idx = search_indices[best_drop_idx + recovery_points[0]]
                fault_end_time = time[fault_end_idx]
                return fault_start_time, fault_end_time
    
    # 3) Final fallback: use code-based calculation
    return code_fault_start, code_fault_end

def plot_all_variables():
    """Plot all variables from the simulation CSV file"""
    
    try:
        # Read the CSV file using robust loader
        df = load_plot_csv()
        time = df['time'].values
        
        print("=== Available Variables ===")
        for i, col in enumerate(df.columns):
            print(f"{i:2d}. {col}")
        
        # Create a comprehensive plot with multiple subplots
        fig, axes = plt.subplots(4, 4, figsize=(20, 16))
        fig.suptitle('Power System Dynamics: Complete Variable Analysis', fontsize=16, fontweight='bold')
        
        # Flatten axes for easier indexing
        axes = axes.flatten()
        
        # Define variable groups with units and descriptions
        variables = [
            ('delta', 'Rotor Angle δ (rad)', 'Generator rotor angle'),
            ('slip', 'Slip s (pu)', 'Generator slip'),
            ('Efd', 'Field Voltage Efd (pu)', 'Exciter field voltage'),
            ('Vref', 'Reference Voltage (pu)', 'AVR reference voltage') if 'Vref' in df.columns else None,
            ('Ed2', 'Ed″ (pu)', 'Sub-transient d-axis voltage'),
            ('Eq2', 'Eq″ (pu)', 'Sub-transient q-axis voltage'), 
            ('Ed1', 'Ed\' (pu)', 'Transient d-axis voltage'),
            ('Eq1', 'Eq\' (pu)', 'Transient q-axis voltage'),
            ('vq', 'VQ (pu)', 'Q-axis terminal voltage'),
            ('vd', 'VD (pu)', 'D-axis terminal voltage'),
            ('E_dummy', 'E_dummy (pu)', 'Dummy voltage'),
            ('Vt', 'Vt (pu)', 'Terminal voltage magnitude'),
            ('iq_0', 'iq (pu)', 'Q-axis current'),
            ('id_0', 'id (pu)', 'D-axis current'),
            ('mech_power', 'Pm (pu)', 'Mechanical power'),
            ('elec_power', 'Pe (pu)', 'Electrical power')
        ]
        
        # Filter out None entries
        variables = [v for v in variables if v]
        
        # Plot each variable
        for i, (var_name, ylabel, description) in enumerate(variables):
            if var_name in df.columns and i < len(axes):
                data = df[var_name].values
                
                axes[i].plot(time, data, linewidth=1.5, color=f'C{i%10}')
                axes[i].set_title(f'{description}\n{var_name}', fontsize=10, fontweight='bold')
                axes[i].set_xlabel('Time (s)')
                axes[i].set_ylabel(ylabel)
                axes[i].grid(True, alpha=0.3)
                
                # Add some statistics as text
                axes[i].text(0.02, 0.98, f'Min: {np.min(data):.4f}\nMax: {np.max(data):.4f}', 
                           transform=axes[i].transAxes, fontsize=8, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Remove empty subplots
        for i in range(len(variables), len(axes)):
            fig.delaxes(axes[i])
        
        plt.tight_layout()
        plt.savefig('report_plot/all_variables_vs_time.png', dpi=300, bbox_inches='tight')
        print("\nComprehensive plot saved as 'report_plot/all_variables_vs_time.png'")
        
        # Create additional focused plots for key variables
        create_key_plots(df, time)
        
    except FileNotFoundError:
        print("Error: PLOT_SAVED_BUS.csv file not found!")
        print("Please run the simulation first using './test'")
    except Exception as e:
        print(f"Error reading/plotting data: {e}")

def create_key_plots(df, time):
    """Create focused plots for key variable groups"""
    
    try:
        print("Creating generator dynamics plot...")
        # 1. Generator Dynamics (Delta and Slip)
        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Delta in degrees
        delta_deg = df['delta'].values * 180 / np.pi
        ax1.plot(time, delta_deg, 'b-', linewidth=2, label='Rotor Angle δ')
        ax1.set_ylabel('Rotor Angle δ (degrees)')
        ax1.set_title('Generator Rotor Dynamics')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Slip
        ax2.plot(time, df['slip'].values, 'r-', linewidth=2, label='Slip s')
        ax2.set_xlabel('Time (seconds)')
        ax2.set_ylabel('Slip s (pu)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig('report_plot/generator_dynamics.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Generator dynamics plot saved as 'report_plot/generator_dynamics.png'")
    except Exception as e:
        print(f"❌ Generator dynamics plot failed: {e}")
    
    try:
        print("Creating voltage analysis plot...")
        # 2. Voltage Analysis
        fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Terminal voltages
        ax1.plot(time, df['Vt'].values, 'g-', linewidth=2, label='Terminal Voltage Vt')
        # Calculate reference voltage estimate (assume constant 1.04 pu for typical system)
        vref_est = np.ones_like(time) * 1.04
        ax1.plot(time, vref_est, 'g--', linewidth=2, label='Reference Voltage Vref (est.)')
        ax1.set_ylabel('Voltage (pu)')
        ax1.set_title('Terminal and Reference Voltages')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # dq voltages
        ax2.plot(time, df['vd'].values, 'b-', linewidth=2, label='VD')
        ax2.plot(time, df['vq'].values, 'r-', linewidth=2, label='VQ')
        ax2.set_ylabel('Voltage (pu)')
        ax2.set_title('dq-axis Terminal Voltages')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # Internal voltages
        ax3.plot(time, df['Eq1'].values, 'c-', linewidth=2, label="Eq'")
        ax3.plot(time, df['Ed1'].values, 'm-', linewidth=2, label="Ed'")
        ax3.set_ylabel('Voltage (pu)')
        ax3.set_title('Transient Internal Voltages')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # Field voltage and error (simplified - no twin axes)
        ax4.plot(time, df['Efd'].values, 'orange', linewidth=2, label='Field Voltage Efd')
        ax4.set_xlabel('Time (seconds)')
        ax4.set_ylabel('Field Voltage (pu)')
        ax4.set_title('Exciter Field Voltage')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        
        plt.tight_layout()
        plt.savefig('report_plot/voltage_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Voltage analysis plot saved as 'report_plot/voltage_analysis.png'")
    except Exception as e:
        print(f"❌ Voltage analysis plot failed: {e}")
    
    try:
        print("Creating power analysis plot...")
        # 3. Power Analysis
        fig3, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # Power comparison
        ax1.plot(time, df['mech_power'].values, 'b-', linewidth=2, label='Mechanical Power Pm')
        ax1.plot(time, df['elec_power'].values, 'r-', linewidth=2, label='Electrical Power Pe')
        ax1.set_ylabel('Power (pu)')
        ax1.set_title('Mechanical vs Electrical Power')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Power difference (acceleration power)
        power_diff = df['mech_power'].values - df['elec_power'].values
        ax2.plot(time, power_diff, 'g-', linewidth=2, label='Pm - Pe (Acceleration Power)')
        ax2.set_xlabel('Time (seconds)')
        ax2.set_ylabel('Power Difference (pu)')
        ax2.set_title('Acceleration Power (Pm - Pe)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig('report_plot/power_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Power analysis plot saved as 'report_plot/power_analysis.png'")
    except Exception as e:
        print(f"❌ Power analysis plot failed: {e}")
    
    try:
        print("Creating current analysis plot...")
        # 4. Current Analysis
        fig4, ax = plt.subplots(1, 1, figsize=(12, 6))
        
        ax.plot(time, df['id_0'].values, 'b-', linewidth=2, label='d-axis Current id')
        ax.plot(time, df['iq_0'].values, 'r-', linewidth=2, label='q-axis Current iq')
        
        # Calculate total current magnitude
        i_total = np.sqrt(df['id_0'].values**2 + df['iq_0'].values**2)
        ax.plot(time, i_total, 'g-', linewidth=2, label='Total Current |I|')
        
        ax.set_xlabel('Time (seconds)')
        ax.set_ylabel('Current (pu)')
        ax.set_title('Generator Current Components')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        plt.savefig('report_plot/current_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Current analysis plot saved as 'report_plot/current_analysis.png'")
    except Exception as e:
        print(f"❌ Current analysis plot failed: {e}")

def plot_fault_period():
    """Plot the fault period showing the dramatic voltage and power changes"""
    
    # Read simulation data using robust loader
    df = load_plot_csv()
    time = df['time'].values
    
    # Detect fault timing from data or code
    fault_start, fault_end = _detect_fault_timing_from_data(df, time)
    
    # Create time window around fault (±10 seconds)
    window_start = fault_start - 10
    window_end = fault_end + 10
    
    # Find indices for the window
    start_idx = np.argmin(np.abs(time - window_start))
    end_idx = np.argmin(np.abs(time - window_end))
    
    # Extract windowed data
    time_window = time[start_idx:end_idx]
    vt_window = df['Vt'].values[start_idx:end_idx]
    pe_window = df['elec_power'].values[start_idx:end_idx]
    delta_window = df['delta'].values[start_idx:end_idx] * 180 / np.pi  # Convert to degrees
    
    # Create the plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10))
    
    # Plot 1: Terminal Voltage
    ax1.plot(time_window, vt_window, 'b-', linewidth=2, label='Terminal Voltage Vt')
    ax1.axvspan(fault_start, fault_end, alpha=0.3, color='red', label='FAULT PERIOD')
    ax1.axvline(fault_start, color='red', linestyle='--', alpha=0.7, label='Fault Start')
    ax1.axvline(fault_end, color='green', linestyle='--', alpha=0.7, label='Fault Cleared')
    ax1.set_ylabel('Terminal Voltage (pu)')
    ax1.set_title('POWER SYSTEM FAULT: Terminal Voltage Drop & Recovery')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 1.2)
    
    # Add annotations
    pre_fault_vt = np.mean(vt_window[time_window < fault_start])
    fault_vt = np.mean(vt_window[(time_window >= fault_start) & (time_window <= fault_end)])
    post_fault_vt = np.mean(vt_window[time_window > fault_end])
    
    ax1.text(fault_start - 5, 1.1, f'Pre-fault:\n{pre_fault_vt:.3f} pu', 
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax1.text(fault_start + 0.05, 0.3, f'During fault:\n{fault_vt:.3f} pu\n({(1-fault_vt/pre_fault_vt)*100:.0f}% drop)', 
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    ax1.text(fault_end + 2, 1.1, f'Post-fault:\n{post_fault_vt:.3f} pu', 
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # Plot 2: Electrical Power
    ax2.plot(time_window, pe_window, 'g-', linewidth=2, label='Electrical Power Pe')
    ax2.axvspan(fault_start, fault_end, alpha=0.3, color='red', label='FAULT PERIOD')
    ax2.axvline(fault_start, color='red', linestyle='--', alpha=0.7)
    ax2.axvline(fault_end, color='green', linestyle='--', alpha=0.7)
    ax2.set_ylabel('Electrical Power (pu)')
    ax2.set_title('Electrical Power During Fault')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: Rotor Angle
    ax3.plot(time_window, delta_window, 'r-', linewidth=2, label='Rotor Angle δ')
    ax3.axvspan(fault_start, fault_end, alpha=0.3, color='red', label='FAULT PERIOD')
    ax3.axvline(fault_start, color='red', linestyle='--', alpha=0.7)
    ax3.axvline(fault_end, color='green', linestyle='--', alpha=0.7)
    ax3.set_xlabel('Time (seconds)')
    ax3.set_ylabel('Rotor Angle (degrees)')
    ax3.set_title('Generator Rotor Angle Response')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Highlight the exact fault duration
    fault_duration = fault_end - fault_start
    fault_cycles_actual = fault_duration * 60.0  # Convert seconds to cycles
    fig.suptitle(f'FAULT ANALYSIS: Bus 0 Terminal Fault\n'
                f'Duration: {fault_duration:.3f} seconds ({fault_cycles_actual:.1f} cycles)\n'
                f'Applied: {fault_start:.1f}s → Cleared: {fault_end:.3f}s', 
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('report_plot/fault_period_analysis.png', dpi=300, bbox_inches='tight')
    print("✅ Fault period plot saved as 'report_plot/fault_period_analysis.png'")

def plot_voltage_recovery():
    """Create detailed voltage recovery analysis plot"""
    
    # Read simulation data using robust loader
    df = load_plot_csv()
    time = df['time'].values
    vt = df['Vt'].values
    
    # Detect fault timing from data
    fault_start, fault_end = _detect_fault_timing_from_data(df, time)
    
    # Analysis periods
    pre_fault_start = fault_start - 10.0
    pre_fault_end = fault_start - 0.1
    post_fault_start = fault_end + 5.0
    post_fault_end = fault_end + 100.0
    
    # Find indices
    pre_fault_mask = (time >= pre_fault_start) & (time <= pre_fault_end)
    post_fault_mask = (time >= post_fault_start) & (time <= post_fault_end)
    
    # Calculate averages
    pre_fault_avg = np.mean(vt[pre_fault_mask])
    post_fault_avg = np.mean(vt[post_fault_mask])
    recovery_threshold = 0.95 * pre_fault_avg
    
    # Create detailed plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))
    
    # Plot 1: Full time series
    ax1.plot(time, vt, 'b-', linewidth=1, label='Terminal Voltage Vt')
    ax1.axvspan(fault_start, fault_end, alpha=0.3, color='red', label='FAULT PERIOD')
    ax1.axvline(fault_start, color='red', linestyle='--', alpha=0.7)
    ax1.axvline(fault_end, color='green', linestyle='--', alpha=0.7)
    
    # Add reference lines
    ax1.axhline(pre_fault_avg, color='blue', linestyle=':', alpha=0.7, label=f'Pre-fault avg: {pre_fault_avg:.3f} pu')
    ax1.axhline(recovery_threshold, color='orange', linestyle=':', alpha=0.7, label=f'Recovery threshold: {recovery_threshold:.3f} pu')
    
    ax1.set_xlabel('Time (seconds)')
    ax1.set_ylabel('Terminal Voltage (pu)')
    ax1.set_title('Terminal Voltage Recovery Analysis - Full Simulation')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim(0, 1.2)
    
    # Plot 2: Zoomed view around fault
    fault_window = 20  # ±20 seconds around fault
    window_mask = (time >= fault_start - fault_window) & (time <= fault_end + fault_window)
    
    ax2.plot(time[window_mask], vt[window_mask], 'b-', linewidth=2, label='Terminal Voltage Vt')
    ax2.axvspan(fault_start, fault_end, alpha=0.3, color='red', label='FAULT PERIOD')
    ax2.axvline(fault_start, color='red', linestyle='--', alpha=0.7, label='Fault Start')
    ax2.axvline(fault_end, color='green', linestyle='--', alpha=0.7, label='Fault Cleared')
    
    ax2.axhline(pre_fault_avg, color='blue', linestyle=':', alpha=0.7, label=f'Pre-fault: {pre_fault_avg:.3f} pu')
    
    ax2.set_xlabel('Time (seconds)')
    ax2.set_ylabel('Terminal Voltage (pu)')
    ax2.set_title('Terminal Voltage Recovery - Detailed View Around Fault')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('report_plot/voltage_recovery_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Detailed voltage recovery plot saved as 'report_plot/voltage_recovery_analysis.png'")

def plot_stability_overview():
    """Create stability overview plots for comprehensive assessment"""
    
    # Read simulation data using robust loader
    df = load_plot_csv()
    time = df['time'].values
    
    # Detect fault timing from data
    fault_start, fault_end = _detect_fault_timing_from_data(df, time)
    
    # Create comprehensive stability plots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Full time series of key variables
    ax1 = plt.subplot(3, 2, 1)
    plt.plot(time, df['Vt'], 'b-', linewidth=1, label='Terminal Voltage')
    plt.axvspan(fault_start, fault_end, alpha=0.3, color='red')
    plt.ylabel('Voltage (pu)')
    plt.title('Terminal Voltage - Full Simulation')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    ax2 = plt.subplot(3, 2, 2)
    plt.plot(time, df['delta'] * 180 / np.pi, 'r-', linewidth=1, label='Rotor Angle')
    plt.axvspan(fault_start, fault_end, alpha=0.3, color='red')
    plt.ylabel('Angle (degrees)')
    plt.title('Rotor Angle - Full Simulation')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    ax3 = plt.subplot(3, 2, 3)
    plt.plot(time, df['slip'], 'g-', linewidth=1, label='Slip')
    plt.axvspan(fault_start, fault_end, alpha=0.3, color='red')
    plt.ylabel('Slip (pu)')
    plt.title('Generator Slip - Full Simulation')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    ax4 = plt.subplot(3, 2, 4)
    plt.plot(time, df['elec_power'], 'purple', linewidth=1, label='Electrical Power')
    plt.plot(time, df['mech_power'], 'orange', linewidth=1, label='Mechanical Power')
    plt.axvspan(fault_start, fault_end, alpha=0.3, color='red')
    plt.ylabel('Power (pu)')
    plt.title('Power Balance')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Plot 5: Pre-fault stability (first 60 seconds)
    ax5 = plt.subplot(3, 2, 5)
    mask_60s = time <= 60
    plt.plot(time[mask_60s], df['Vt'][mask_60s], 'b-', linewidth=2)
    plt.ylabel('Voltage (pu)')
    plt.xlabel('Time (s)')
    plt.title('Initial 60s - Voltage Settling')
    plt.grid(True, alpha=0.3)
    
    # Plot 6: Post-fault stability (last 60 seconds)
    ax6 = plt.subplot(3, 2, 6)
    mask_last60s = time >= (time[-1] - 60)
    plt.plot(time[mask_last60s], df['Vt'][mask_last60s], 'b-', linewidth=2)
    plt.ylabel('Voltage (pu)')
    plt.xlabel('Time (s)')
    plt.title('Final 60s - Final Settling')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('report_plot/stability_overview.png', dpi=300, bbox_inches='tight')
    print("Stability overview plots saved as 'report_plot/stability_overview.png'")

def main():
    """Run all plotting functions"""
    print("="*70)
    print("POWER SYSTEM DYNAMICS - COMPREHENSIVE PLOTTING SUITE")
    print("="*70)
    
    try:
        plot_all_variables()
        plot_fault_period()
        plot_voltage_recovery()
        plot_stability_overview()
        
        print(f"\n{'✅ ALL PLOTS GENERATED SUCCESSFULLY'}")
        print("="*70)
        print("Generated files in report_plot/:")
        print("• all_variables_vs_time.png - All variables overview")
        print("• generator_dynamics.png - Rotor angle and slip")
        print("• voltage_analysis.png - Voltage components")
        print("• power_analysis.png - Power balance")
        print("• current_analysis.png - Current components")
        print("• fault_period_analysis.png - Fault response")
        print("• voltage_recovery_analysis.png - Recovery assessment")
        print("• stability_overview.png - Stability overview")
        print("="*70)
        
    except Exception as e:
        print(f"Error in plotting: {e}")

if __name__ == "__main__":
    main() 