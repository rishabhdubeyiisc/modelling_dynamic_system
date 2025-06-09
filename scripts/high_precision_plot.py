#!/usr/bin/env python3
"""
High Precision Plotting Script for Power System Dynamics
Generates publication-quality plots with enhanced detail and precision
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_high_precision_plots():
    """Generate high-precision plots for all key variables"""
    
    # Read simulation data
    df = pd.read_csv('PLOT_SAVED_BUS.csv', sep=',', header=0)
    
    # Set up high-quality plot parameters
    plt.rcParams.update({
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'lines.linewidth': 1.5,
        'grid.alpha': 0.3
    })
    
    # Create output directory
    os.makedirs('report_plot', exist_ok=True)
    
    time = df['time'].values
    
    # 1. GENERATOR DYNAMICS - High precision
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Rotor angle in degrees
    delta_deg = df['delta'].values * 180 / np.pi
    ax1.plot(time, delta_deg, 'b-', linewidth=1.5, label='Rotor Angle δ')
    ax1.set_ylabel('Rotor Angle δ (degrees)')
    ax1.set_title('Rotor Angle Dynamics (High Precision)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Slip
    ax2.plot(time, df['slip'].values * 1000, 'r-', linewidth=1.5, label='Slip (×1000)')
    ax2.set_ylabel('Slip (×1000 pu)')
    ax2.set_title('Generator Slip (High Precision)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Field voltage
    ax3.plot(time, df['Efd'].values, 'g-', linewidth=1.5, label='Field Voltage Efd')
    ax3.set_ylabel('Field Voltage Efd (pu)')
    ax3.set_xlabel('Time (s)')
    ax3.set_title('Exciter Field Voltage')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Terminal voltage
    ax4.plot(time, df['Vt'].values, 'm-', linewidth=1.5, label='Terminal Voltage Vt')
    ax4.set_ylabel('Terminal Voltage Vt (pu)')
    ax4.set_xlabel('Time (s)')
    ax4.set_title('Terminal Voltage')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('report_plot/high_precision_generator_dynamics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. EMF DYNAMICS - All EMF components
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    ax1.plot(time, df['Eq1'].values, 'b-', linewidth=1.5, label="Eq' (transient)")
    ax1.plot(time, df['Eq2'].values, 'r-', linewidth=1.5, label='Eq" (sub-transient)')
    ax1.set_ylabel('Q-axis EMF (pu)')
    ax1.set_title('Q-axis EMF Components')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.plot(time, df['Ed1'].values, 'b-', linewidth=1.5, label="Ed' (transient)")
    ax2.plot(time, df['Ed2'].values, 'r-', linewidth=1.5, label='Ed" (sub-transient)')
    ax2.set_ylabel('D-axis EMF (pu)')
    ax2.set_title('D-axis EMF Components')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # EMF magnitude
    eq_mag = np.sqrt(df['Eq1'].values**2 + df['Ed1'].values**2)
    ax3.plot(time, eq_mag, 'g-', linewidth=1.5, label="EMF Magnitude |E'|")
    ax3.set_ylabel('EMF Magnitude (pu)')
    ax3.set_xlabel('Time (s)')
    ax3.set_title('Transient EMF Magnitude')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # E_dummy (dummy coil EMF)
    ax4.plot(time, df['E_dummy'].values, 'orange', linewidth=1.5, label='E_dummy')
    ax4.set_ylabel('E_dummy (pu)')
    ax4.set_xlabel('Time (s)')
    ax4.set_title('Dummy Coil EMF')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('report_plot/high_precision_emf_dynamics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. CURRENT DYNAMICS - Detailed current analysis
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # D-Q currents
    ax1.plot(time, df['id_0'].values, 'b-', linewidth=1.5, label='D-axis current id')
    ax1.plot(time, df['iq_0'].values, 'r-', linewidth=1.5, label='Q-axis current iq')
    ax1.set_ylabel('Current (pu)')
    ax1.set_title('D-Q Axis Currents')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Current magnitude
    i_mag = np.sqrt(df['id_0'].values**2 + df['iq_0'].values**2)
    ax2.plot(time, i_mag, 'g-', linewidth=1.5, label='Current Magnitude |I|')
    ax2.set_ylabel('Current Magnitude (pu)')
    ax2.set_title('Total Current Magnitude')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # D-Q voltages at terminals
    ax3.plot(time, df['vd'].values, 'b-', linewidth=1.5, label='D-axis voltage vd')
    ax3.plot(time, df['vq'].values, 'r-', linewidth=1.5, label='Q-axis voltage vq')
    ax3.set_ylabel('Voltage (pu)')
    ax3.set_xlabel('Time (s)')
    ax3.set_title('D-Q Terminal Voltages')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Power factor and current angle
    current_angle = np.arctan2(df['iq_0'].values, df['id_0'].values) * 180 / np.pi
    ax4.plot(time, current_angle, 'purple', linewidth=1.5, label='Current Angle')
    ax4.set_ylabel('Current Angle (degrees)')
    ax4.set_xlabel('Time (s)')
    ax4.set_title('Current Phase Angle')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('report_plot/high_precision_current_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. POWER ANALYSIS - Detailed power dynamics
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Mechanical vs Electrical power
    ax1.plot(time, df['mech_power'].values, 'b-', linewidth=1.5, label='Mechanical Power Pm')
    ax1.plot(time, df['elec_power'].values, 'r-', linewidth=1.5, label='Electrical Power Pe')
    ax1.set_ylabel('Power (pu)')
    ax1.set_title('Mechanical vs Electrical Power')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Power imbalance (high precision)
    power_imbalance = df['mech_power'].values - df['elec_power'].values
    ax2.plot(time, power_imbalance * 1000, 'g-', linewidth=1.5, label='Power Imbalance (×1000)')
    ax2.set_ylabel('Power Imbalance (×1000 pu)')
    ax2.set_title('Power Imbalance (High Precision)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Power vs speed correlation
    ax3.scatter(df['slip'].values * 1000, power_imbalance * 1000, alpha=0.6, s=1)
    ax3.set_xlabel('Slip (×1000 pu)')
    ax3.set_ylabel('Power Imbalance (×1000 pu)')
    ax3.set_title('Power Imbalance vs Slip')
    ax3.grid(True, alpha=0.3)
    
    # Running average power
    window = min(100, len(time)//10)
    if window > 1:
        pm_avg = pd.Series(df['mech_power'].values).rolling(window).mean()
        pe_avg = pd.Series(df['elec_power'].values).rolling(window).mean()
        ax4.plot(time, pm_avg, 'b-', linewidth=1.5, label=f'Pm (avg {window} pts)')
        ax4.plot(time, pe_avg, 'r-', linewidth=1.5, label=f'Pe (avg {window} pts)')
        ax4.set_ylabel('Power (pu)')
        ax4.set_xlabel('Time (s)')
        ax4.set_title('Smoothed Power Curves')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
    
    plt.tight_layout()
    plt.savefig('report_plot/high_precision_power_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. ZOOM PLOTS - Critical periods
    # Find interesting periods (e.g., when slip changes significantly)
    slip_diff = np.abs(np.diff(df['slip'].values))
    if len(slip_diff) > 0:
        max_change_idx = np.argmax(slip_diff)
        critical_time = time[max_change_idx]
        
        # Plot 10 seconds around critical event
        window_start = max(0, critical_time - 5)
        window_end = min(time[-1], critical_time + 5)
        
        mask = (time >= window_start) & (time <= window_end)
        time_zoom = time[mask]
        
        if len(time_zoom) > 10:  # Only if we have enough data points
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            
            # Zoomed rotor dynamics
            ax1.plot(time_zoom, df['delta'].values[mask] * 180 / np.pi, 'b-', linewidth=2)
            ax1.set_ylabel('Rotor Angle (deg)')
            ax1.set_title(f'Rotor Angle (Zoom: {window_start:.1f}-{window_end:.1f}s)')
            ax1.grid(True, alpha=0.3)
            
            # Zoomed slip
            ax2.plot(time_zoom, df['slip'].values[mask] * 1000, 'r-', linewidth=2)
            ax2.set_ylabel('Slip (×1000 pu)')
            ax2.set_title(f'Slip (Zoom: {window_start:.1f}-{window_end:.1f}s)')
            ax2.grid(True, alpha=0.3)
            
            # Zoomed voltage
            ax3.plot(time_zoom, df['Vt'].values[mask], 'g-', linewidth=2)
            ax3.set_ylabel('Terminal Voltage (pu)')
            ax3.set_xlabel('Time (s)')
            ax3.set_title(f'Terminal Voltage (Zoom: {window_start:.1f}-{window_end:.1f}s)')
            ax3.grid(True, alpha=0.3)
            
            # Zoomed power
            ax4.plot(time_zoom, df['mech_power'].values[mask], 'b-', linewidth=2, label='Pm')
            ax4.plot(time_zoom, df['elec_power'].values[mask], 'r-', linewidth=2, label='Pe')
            ax4.set_ylabel('Power (pu)')
            ax4.set_xlabel('Time (s)')
            ax4.set_title(f'Power (Zoom: {window_start:.1f}-{window_end:.1f}s)')
            ax4.grid(True, alpha=0.3)
            ax4.legend()
            
            plt.tight_layout()
            plt.savefig('report_plot/high_precision_critical_period.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    print("✅ High-precision plots generated:")
    print("   • high_precision_generator_dynamics.png")
    print("   • high_precision_emf_dynamics.png") 
    print("   • high_precision_current_analysis.png")
    print("   • high_precision_power_analysis.png")
    print("   • high_precision_critical_period.png")

if __name__ == "__main__":
    create_high_precision_plots() 