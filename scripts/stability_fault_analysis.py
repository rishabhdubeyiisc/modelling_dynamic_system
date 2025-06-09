#!/usr/bin/env python3
"""
Power System Dynamics - Comprehensive Analysis Suite
Combines stability analysis, voltage recovery, and fault timing functionality
"""

import pandas as pd
import numpy as np
import re
from scipy import signal

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

def _detect_fault_timing_from_data(df, time):
    """Detect actual fault timing from simulation data by finding voltage drops."""
    # First try code-based calculation (more reliable)
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
    
    # Fallback: use code-based calculation (most reliable)
    print(f"⚠️  Using code-based fault timing (data detection failed)")
    return code_fault_start, code_fault_end

def _detect_reference_voltage(df):
    """Detect reference voltage from simulation data."""
    # Method 1: Try to read Vref column if it exists
    if 'Vref' in df.columns:
        vref_values = df['Vref'].values
        # Use the most common value (mode) as reference
        vref = np.median(vref_values)  # Use median as it's more robust
        return vref
    
    # Method 2: Estimate Vref based on typical power system values
    # For generator buses, Vref is typically around 1.0-1.05 pu
    time = df['time'].values
    vt = df['Vt'].values
    
    # Use average of first 10% of simulation, but bound it to realistic values
    pre_fault_end_idx = int(0.1 * len(time))
    vt_initial_avg = np.mean(vt[:pre_fault_end_idx])
    
    # For power system stability, Vref is typically close to rated voltage (1.0 pu)
    # Bound the estimate to reasonable values (0.95 to 1.1 pu)
    vref = max(0.95, min(1.1, vt_initial_avg))
    
    return vref

def check_voltage_recovery():
    """Check terminal voltage recovery after fault clearing"""
    
    # Read simulation data
    df = pd.read_csv('PLOT_SAVED_BUS.csv', sep=',', header=0)
    time = df['time'].values
    vt = df['Vt'].values
    
    # Detect fault timing from data
    fault_start, fault_end = _detect_fault_timing_from_data(df, time)
    print(f"🔍 Detected fault timing: {fault_start:.3f}s → {fault_end:.3f}s")
    
    # Adaptive analysis periods based on simulation duration
    sim_duration = time[-1] - time[0]
    
    # Pre-fault: 10 seconds before fault (or 5% of pre-fault time, whichever is smaller)
    pre_fault_window = min(10.0, fault_start * 0.05)
    pre_fault_start = fault_start - pre_fault_window
    pre_fault_end = fault_start - 0.1
    
    # Post-fault: adaptive window based on remaining simulation time  
    remaining_time = sim_duration - fault_end
    post_fault_start = fault_end + min(5.0, remaining_time * 0.05)
    post_fault_end = fault_end + min(100.0, remaining_time * 0.5)
    
    # Find indices
    pre_fault_mask = (time >= pre_fault_start) & (time <= pre_fault_end)
    fault_mask = (time >= fault_start) & (time <= fault_end)
    post_fault_mask = (time >= post_fault_start) & (time <= post_fault_end)
    
    # Calculate averages
    pre_fault_avg = np.mean(vt[pre_fault_mask])
    fault_avg = np.mean(vt[fault_mask])
    post_fault_avg = np.mean(vt[post_fault_mask])
    
    # Final recovery check (last 10% of simulation or 50 seconds, whichever is smaller)
    final_window = min(50.0, sim_duration * 0.1)
    final_start = time[-1] - final_window
    final_period_mask = (time >= final_start) & (time <= time[-1])
    final_avg = np.mean(vt[final_period_mask])
    final_std = np.std(vt[final_period_mask])
    
    print("=" * 70)
    print("TERMINAL VOLTAGE RECOVERY ANALYSIS")
    print("=" * 70)
    print(f"Pre-fault voltage (t={pre_fault_start}-{pre_fault_end}s):     {pre_fault_avg:.4f} pu")
    print(f"During fault voltage (t={fault_start}-{fault_end}s):         {fault_avg:.4f} pu")
    print(f"Post-fault voltage (t={post_fault_start}-{post_fault_end}s):  {post_fault_avg:.4f} pu")
    print(f"Final voltage (t={final_start:.1f}-{time[-1]:.1f}s):                 {final_avg:.4f} ± {final_std:.6f} pu")
    print()
    
    # Recovery metrics
    voltage_drop = pre_fault_avg - fault_avg
    immediate_recovery = post_fault_avg - fault_avg
    full_recovery = final_avg - fault_avg
    recovery_percentage = (final_avg / pre_fault_avg) * 100
    
    print("RECOVERY METRICS:")
    print(f"Voltage drop during fault:       {voltage_drop:.4f} pu ({(voltage_drop/pre_fault_avg)*100:.1f}%)")
    print(f"Immediate recovery (5s after):  {immediate_recovery:.4f} pu")
    print(f"Full recovery:                   {full_recovery:.4f} pu")
    print(f"Final voltage vs pre-fault:      {recovery_percentage:.1f}%")
    print()
    
    # Check if voltage has recovered
    recovery_threshold = 0.95 * pre_fault_avg  # 95% recovery considered good
    
    if final_avg >= recovery_threshold and final_std < 0.001:
        print("✅ VOLTAGE RECOVERY: SUCCESS")
        print(f"   Terminal voltage has recovered to {recovery_percentage:.1f}% of pre-fault value")
        print(f"   System is stable (std dev = {final_std:.6f})")
        recovery_status = "SUCCESS"
    elif final_avg >= recovery_threshold:
        print("⚠️  VOLTAGE RECOVERY: PARTIAL")
        print(f"   Voltage level recovered but still oscillating (std dev = {final_std:.6f})")
        recovery_status = "PARTIAL"
    else:
        print("❌ VOLTAGE RECOVERY: FAILED")
        print(f"   Final voltage only {recovery_percentage:.1f}% of pre-fault value")
        recovery_status = "FAILED"
    
    print("=" * 70)
    
    return recovery_status, {
        'pre_fault': pre_fault_avg,
        'fault': fault_avg,
        'post_fault': post_fault_avg,
        'final': final_avg,
        'recovery_percent': recovery_percentage
    }

def analyze_fault_timing():
    """Analyze fault timing based on code parameters"""
    
    print("\n📋 FAULT TIMING ANALYSIS:")
    print("-" * 40)
    
    # Parameters from the code
    FAULT_CYCLES, CLEAR_DISTURBANCES_SEC = _read_fault_parameters_from_main()
    PERIOD = 1.0/60.0  # seconds (pow(60,-1))
    FAULTED_TIME = FAULT_CYCLES * PERIOD  # seconds
    
    print(f"FAULT_CYCLES:           {FAULT_CYCLES} cycles")
    print(f"CLEAR_DISTURBANCES:     {CLEAR_DISTURBANCES_SEC} seconds ({CLEAR_DISTURBANCES_SEC/60:.1f} minutes)")
    print(f"PERIOD:                 {PERIOD:.6f} seconds (1/60 Hz)")
    print(f"FAULTED_TIME:           {FAULTED_TIME:.6f} seconds")
    
    fault_start_time = 2 * CLEAR_DISTURBANCES_SEC
    fault_end_time = fault_start_time + FAULTED_TIME
    
    print(f"Fault starts at:        t > {fault_start_time:.3f} seconds ({fault_start_time/60:.1f} minutes)")
    print(f"Fault active during:    {fault_start_time:.3f} < t < {fault_end_time:.3f} seconds")
    print(f"Fault ends at:          t > {fault_end_time:.3f} seconds")
    print(f"Fault duration:         {FAULTED_TIME:.6f} seconds = {FAULT_CYCLES} cycles")
    
    return fault_start_time, fault_end_time

def analyze_stability():
    """Comprehensive stability analysis of the power system"""
    
    # Read simulation data
    df = pd.read_csv('PLOT_SAVED_BUS.csv', sep=',', header=0)
    time = df['time'].values
    
    # Detect fault timing and reference voltage
    fault_start, fault_end = _detect_fault_timing_from_data(df, time)
    vref = _detect_reference_voltage(df)
    
    print("\n🔍 STABILITY ANALYSIS:")
    print("-" * 50)
    
    # 1. INITIAL CONDITIONS ANALYSIS
    print("\n📊 Initial Condition Assessment:")
    vt_0 = df['Vt'].iloc[0]
    initial_voltage_error = abs(vref - vt_0)
    print(f"   Terminal Voltage Vt(0):  {vt_0:.6f} pu")
    print(f"   Reference Voltage Vref:  {vref:.6f} pu (detected)")
    print(f"   Initial Voltage Error:   {initial_voltage_error:.6f} pu ({initial_voltage_error/vref*100:.2f}%)")
    
    if initial_voltage_error > 0.01:
        print("   ❌ LARGE INITIAL VOLTAGE ERROR - This will cause instability!")
    else:
        print("   ✅ Initial voltage error is acceptable")
    
    # 2. PRE-FAULT STABILITY ANALYSIS
    pre_fault_mask = time < fault_start
    analyze_period(df[pre_fault_mask], "PRE-FAULT")
    
    # 3. POST-FAULT STABILITY ANALYSIS
    post_fault_mask = time > fault_end
    analyze_period(df[post_fault_mask], "POST-FAULT")
    
    # 4. OSCILLATION ANALYSIS
    analyze_oscillations(df, time, fault_start, fault_end)
    
    # 5. STABILITY RECOMMENDATIONS
    provide_stability_recommendations(df, time, fault_start, fault_end)

def analyze_period(df_period, period_name):
    """Analyze stability for a specific time period"""
    
    if len(df_period) == 0:
        print(f"   No data for {period_name} period")
        return
    
    print(f"\n📈 {period_name} PERIOD ANALYSIS:")
    
    # Calculate statistics for key variables
    variables = ['Vt', 'delta', 'slip', 'elec_power', 'mech_power']
    
    for var in variables:
        if var in df_period.columns:
            data = df_period[var].values
            if var == 'delta':
                data = data * 180 / np.pi  # Convert to degrees
            
            mean_val = np.mean(data)
            std_val = np.std(data)
            max_val = np.max(data)
            min_val = np.min(data)
            drift = data[-1] - data[0]  # Final - Initial
            
            # Check for stability
            is_stable = std_val < 0.01 * abs(mean_val) and abs(drift) < 0.01 * abs(mean_val)
            stability_status = "✅ STABLE" if is_stable else "❌ UNSTABLE"
            
            units = "°" if var == 'delta' else "pu"
            print(f"   {var:12}: Mean={mean_val:8.4f} {units}, "
                  f"Std={std_val:8.4f}, Range=[{min_val:7.4f}, {max_val:7.4f}], "
                  f"Drift={drift:8.4f} {stability_status}")

def analyze_oscillations(df, time, fault_start, fault_end):
    """Analyze oscillatory behavior in the system"""
    
    print(f"\n🔊 OSCILLATION ANALYSIS:")
    
    # Focus on rotor angle for oscillation analysis
    delta = df['delta'].values * 180 / np.pi  # Convert to degrees
    
    # Pre-fault oscillations
    pre_fault_mask = time < fault_start
    if np.sum(pre_fault_mask) > 100:  # Need sufficient data points
        delta_pre = delta[pre_fault_mask]
        time_pre = time[pre_fault_mask]
        
        # Detect oscillations using FFT
        dt = np.mean(np.diff(time_pre))
        freqs, psd = signal.welch(delta_pre, fs=1/dt, nperseg=min(1024, len(delta_pre)//4))
        
        # Find dominant frequencies (excluding DC component)
        dominant_freq_idx = np.argmax(psd[1:]) + 1
        dominant_freq = freqs[dominant_freq_idx]
        
        print(f"   PRE-FAULT Oscillations:")
        print(f"     Rotor angle swing:     {np.max(delta_pre) - np.min(delta_pre):.2f}°")
        print(f"     Dominant frequency:    {dominant_freq:.4f} Hz ({1/dominant_freq:.2f}s period)")
        print(f"     Standard deviation:    {np.std(delta_pre):.4f}°")
    
    # Post-fault oscillations
    post_fault_mask = time > fault_end
    if np.sum(post_fault_mask) > 100:
        delta_post = delta[post_fault_mask]
        time_post = time[post_fault_mask]
        
        dt = np.mean(np.diff(time_post))
        freqs, psd = signal.welch(delta_post, fs=1/dt, nperseg=min(1024, len(delta_post)//4))
        
        dominant_freq_idx = np.argmax(psd[1:]) + 1
        dominant_freq = freqs[dominant_freq_idx]
        
        print(f"   POST-FAULT Oscillations:")
        print(f"     Rotor angle swing:     {np.max(delta_post) - np.min(delta_post):.2f}°")
        print(f"     Dominant frequency:    {dominant_freq:.4f} Hz ({1/dominant_freq:.2f}s period)")
        print(f"     Standard deviation:    {np.std(delta_post):.4f}°")

def provide_stability_recommendations(df, time, fault_start, fault_end):
    """Provide specific recommendations to improve stability"""
    
    print(f"\n💡 STABILITY RECOMMENDATIONS:")
    
    vt_initial = df['Vt'].iloc[0]
    vt_final = df['Vt'].iloc[-1]
    delta_swing = (np.max(df['delta']) - np.min(df['delta'])) * 180 / np.pi
    
    print("   IDENTIFIED ISSUES:")
    
    # Check initial voltage error  
    vref = _detect_reference_voltage(df)
    if abs(vref - vt_initial) > 0.01:
        print(f"   ❌ Large initial voltage error: {abs(vref - vt_initial):.4f} pu")
        print("      → FIX: Adjust initial conditions in data file")
    
    # Check excessive rotor angle swing
    if delta_swing > 100:
        print(f"   ❌ Excessive rotor angle swing: {delta_swing:.1f}°")
        print("      → FIX: Reduce exciter gains (ka) or increase damping")
    
    # Check voltage drift
    voltage_drift = abs(vt_final - vt_initial)
    if voltage_drift > 0.05:
        print(f"   ❌ Voltage drift: {voltage_drift:.4f} pu")
        print("      → FIX: Check exciter parameters and initial conditions")
    
    # Check slip settling
    final_slip = df['slip'].iloc[-1]
    if abs(final_slip) > 0.001:
        print(f"   ❌ Slip not settled: {final_slip:.6f} pu")
        print("      → FIX: System may need more simulation time or damping")
    
    print("\n   RECOMMENDED ACTIONS:")
    print(f"   1. Fix initial conditions: Set Vt(0) = Vref = {vref:.3f} pu")
    print("   2. Reduce exciter gain: Try ka = 1.0 instead of current value")
    print("   3. Increase damping: Add mechanical damping coefficient")
    print("   4. Reduce time step: Try dt = 0.005s for better numerical stability")
    print("   5. Check load flow: Ensure proper steady-state initialization")

def print_detailed_statistics(df, time):
    """Print detailed statistics for all variables"""
    
    print(f"\n📊 DETAILED SIMULATION STATISTICS:")
    print("-" * 50)
    
    print(f"Simulation Duration: {time[0]:.3f} to {time[-1]:.3f} seconds")
    print(f"Time Steps: {len(time)}")
    print(f"Time Step Size: {np.mean(np.diff(time)):.6f} seconds")
    
    # Key variables analysis
    key_vars = ['delta', 'slip', 'Vt', 'Efd', 'mech_power', 'elec_power']
    
    print(f"\n{'Variable':<12} {'Min':<10} {'Max':<10} {'Mean':<10} {'Std':<10} {'Final':<10}")
    print("-" * 65)
    
    for var in key_vars:
        if var in df.columns:
            data = df[var].values
            if var == 'delta':
                # Convert delta to degrees for display
                data = data * 180 / np.pi
                unit = "°"
            else:
                unit = "pu"
                
            print(f"{var:<12} {np.min(data):<10.4f} {np.max(data):<10.4f} "
                  f"{np.mean(data):<10.4f} {np.std(data):<10.4f} {data[-1]:<10.4f}")
    
    # Stability analysis
    print(f"\n🎯 SYSTEM STABILITY SUMMARY:")
    print("-" * 40)
    
    delta_deg = df['delta'].values * 180 / np.pi
    slip = df['slip'].values
    
    print(f"Rotor Angle Swing: {np.max(delta_deg) - np.min(delta_deg):.3f}°")
    print(f"Maximum Slip: {np.max(np.abs(slip)):.6f} pu")
    print(f"Final Slip: {slip[-1]:.6f} pu")
    
    # Check if system is settling
    final_1000_points = min(1000, len(slip))
    final_slip_std = np.std(slip[-final_1000_points:])
    
    if final_slip_std < 1e-6:
        print("System Status: ✅ STABLE (slip settling to zero)")
    elif final_slip_std < 1e-4:
        print("System Status: ⚠️ MARGINALLY STABLE")
    else:
        print("System Status: ❌ UNSTABLE (slip not settling)")
    
    # Power balance
    power_diff = df['mech_power'].values - df['elec_power'].values
    print(f"Average Power Imbalance: {np.mean(np.abs(power_diff)):.6f} pu")

def generate_comprehensive_report(df, time, recovery_status, recovery_metrics, fault_start, fault_end):
    """Generate a comprehensive analysis report and save to file"""
    
    # Detect parameters
    vref = _detect_reference_voltage(df)
    fault_cycles, clear_disturbances = _read_fault_parameters_from_main()
    
    # Calculate key metrics
    rotor_swing = (np.max(df['delta']) - np.min(df['delta'])) * 180 / np.pi
    max_slip = np.max(np.abs(df['slip']))
    voltage_drop = recovery_metrics['pre_fault'] - recovery_metrics['fault']
    recovery_percentage = recovery_metrics['recovery_percent']
    
    # System stability assessment
    delta_deg = df['delta'].values * 180 / np.pi
    slip = df['slip'].values
    final_1000_points = min(1000, len(slip))
    final_slip_std = np.std(slip[-final_1000_points:])
    
    if final_slip_std < 1e-6:
        stability_status = "STABLE"
    elif final_slip_std < 1e-4:
        stability_status = "MARGINALLY STABLE"
    else:
        stability_status = "UNSTABLE"
    
    # Generate report content
    report_content = f"""
{'='*80}
POWER SYSTEM DYNAMICS - COMPREHENSIVE ANALYSIS REPORT
{'='*80}
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

SIMULATION PARAMETERS:
{'='*40}
• Simulation Duration:      {time[0]:.3f} to {time[-1]:.3f} seconds
• Time Steps:              {len(time)}
• Average Time Step:       {np.mean(np.diff(time)):.6f} seconds
• Reference Voltage:       {vref:.6f} pu (detected)
• Fault Cycles (Code):     {fault_cycles} cycles
• Clear Disturbances:      {clear_disturbances} seconds

FAULT ANALYSIS:
{'='*40}
• Fault Start Time:        {fault_start:.3f} seconds
• Fault End Time:          {fault_end:.3f} seconds  
• Fault Duration:          {fault_end-fault_start:.6f} seconds
• Fault Cycles (Actual):   {(fault_end-fault_start)*60:.1f} cycles

VOLTAGE RECOVERY ANALYSIS:
{'='*40}
• Pre-fault Voltage:       {recovery_metrics['pre_fault']:.4f} pu
• During Fault Voltage:    {recovery_metrics['fault']:.4f} pu
• Post-fault Voltage:      {recovery_metrics['post_fault']:.4f} pu
• Final Voltage:           {recovery_metrics['final']:.4f} pu
• Voltage Drop:            {voltage_drop:.4f} pu ({(voltage_drop/recovery_metrics['pre_fault'])*100:.1f}%)
• Recovery Percentage:     {recovery_percentage:.1f}%
• Recovery Status:         {recovery_status}

ROTOR DYNAMICS:
{'='*40}
• Rotor Angle Swing:       {rotor_swing:.3f}°
• Maximum Slip:            {max_slip:.6f} pu
• Final Slip:              {slip[-1]:.6f} pu
• Final Slip Std Dev:      {final_slip_std:.2e} pu

SYSTEM STABILITY ASSESSMENT:
{'='*40}
• Overall Status:          {stability_status}
• Voltage Recovery:        {'✓' if recovery_status == 'SUCCESS' else '✗'}
• Rotor Stability:         {'✓' if rotor_swing < 100 else '✗'}
• Slip Settling:           {'✓' if abs(slip[-1]) < 0.001 else '✗'}

VARIABLE STATISTICS:
{'='*40}
"""
    
    # Add variable statistics
    key_vars = ['delta', 'slip', 'Vt', 'Efd', 'mech_power', 'elec_power']
    
    report_content += f"{'Variable':<12} {'Min':<10} {'Max':<10} {'Mean':<10} {'Std':<10} {'Final':<10}\n"
    report_content += "-" * 65 + "\n"
    
    for var in key_vars:
        if var in df.columns:
            data = df[var].values
            if var == 'delta':
                data = data * 180 / np.pi
                unit = "°"
            else:
                unit = "pu"
                
            report_content += f"{var:<12} {np.min(data):<10.4f} {np.max(data):<10.4f} "
            report_content += f"{np.mean(data):<10.4f} {np.std(data):<10.4f} {data[-1]:<10.4f}\n"
    
    # Add recommendations
    vt_initial = df['Vt'].iloc[0]
    vt_final = df['Vt'].iloc[-1]
    initial_voltage_error = abs(vref - vt_initial)
    voltage_drift = abs(vt_final - vt_initial)
    
    report_content += f"""
RECOMMENDATIONS:
{'='*40}
"""
    
    if initial_voltage_error > 0.01:
        report_content += f"• CRITICAL: Large initial voltage error ({initial_voltage_error:.4f} pu)\n"
        report_content += f"  → Fix: Set Vt(0) = Vref = {vref:.3f} pu in initial conditions\n"
    
    if rotor_swing > 100:
        report_content += f"• WARNING: Excessive rotor angle swing ({rotor_swing:.1f}°)\n"
        report_content += f"  → Fix: Reduce exciter gains (ka) or increase damping\n"
    
    if voltage_drift > 0.05:
        report_content += f"• WARNING: Voltage drift ({voltage_drift:.4f} pu)\n"
        report_content += f"  → Fix: Check exciter parameters and initial conditions\n"
    
    if abs(slip[-1]) > 0.001:
        report_content += f"• WARNING: Slip not settled ({slip[-1]:.6f} pu)\n"
        report_content += f"  → Fix: System may need more simulation time or damping\n"
    
    if recovery_status != "SUCCESS":
        report_content += f"• CRITICAL: Voltage recovery {recovery_status.lower()}\n"
        report_content += f"  → Fix: Check fault clearing logic and system parameters\n"
    
    if initial_voltage_error <= 0.01 and rotor_swing <= 100 and voltage_drift <= 0.05 and abs(slip[-1]) <= 0.001:
        report_content += f"• ✓ System appears to be operating within acceptable limits\n"
    
    report_content += f"""
FILES GENERATED:
{'='*40}
• stability_analysis_report.txt - This comprehensive report
• All analysis plots saved in report_plot/ folder

ANALYSIS COMPLETE
{'='*40}
Report generated successfully.
"""
    
    # Save report to file
    with open('report_plot/stability_analysis_report.txt', 'w') as f:
        f.write(report_content)
    
    print(f"\n📄 Comprehensive report saved as 'report_plot/stability_analysis_report.txt'")
    return report_content

def main():
    """Run comprehensive analysis suite"""
    print("="*70)
    print("POWER SYSTEM DYNAMICS - COMPREHENSIVE ANALYSIS SUITE")
    print("="*70)
    
    try:
        # Read simulation data first to check if it exists
        df = pd.read_csv('PLOT_SAVED_BUS.csv', sep=',', header=0)
        time = df['time'].values
        
        # Run all analyses
        recovery_status, recovery_metrics = check_voltage_recovery()
        fault_start, fault_end = analyze_fault_timing()
        analyze_stability()
        print_detailed_statistics(df, time)
        
        # Generate comprehensive report
        generate_comprehensive_report(df, time, recovery_status, recovery_metrics, fault_start, fault_end)
        
        print(f"\n{'='*70}")
        print("ANALYSIS COMPLETE")
        print(f"{'='*70}")
        print(f"Voltage Recovery: {recovery_status}")
        print(f"Fault Duration: {fault_end-fault_start:.6f} seconds")
        print(f"Report and analysis saved in report_plot/ folder")
        
    except FileNotFoundError:
        print("❌ Error: PLOT_SAVED_BUS.csv file not found!")
        print("Please run the simulation first using './test'")
    except Exception as e:
        print(f"❌ Error in analysis: {e}")

if __name__ == "__main__":
    main() 