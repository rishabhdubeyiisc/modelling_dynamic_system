#!/usr/bin/env python3
"""
Jacobian Solver Debug Analyzer
Analyzes the detailed per-generator log files to identify convergence issues
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def parse_generator_log(log_file):
    """Parse a generator debug log file and extract convergence data"""
    data = {
        'times': [],
        'iterations': [],
        'residuals': {'Ed': [], 'Eq': [], "Ed''": [], "Eq''": [], 'delta': [], 'slip': [], 'Efd': []},
        'states': {'Ed': [], 'Eq': [], "Ed''": [], "Eq''": [], 'delta': [], 'slip': [], 'Efd': []},
        'currents': {'id': [], 'iq': []},
        'voltages': {'vd': [], 'vq': [], 'Vt': [], 'Vref': [], 'VD': [], 'VQ': []},
        'power': {'Pm': [], 'Pe': []},
        'equation_components': {
            'Ed_time': [], 'Ed_steady': [],
            'Eq_time': [], 'Eq_steady': [],
            'slip_time': [], 'torque_imb': [], 'damping': [],
            'Efd_time': [], 'Efd_steady': []
        },
        'max_residuals': []
    }
    
    if not os.path.exists(log_file):
        print(f"Warning: Log file {log_file} not found")
        return data
    
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Parse header line: --- t=0.000000 iter=0 Gen_0 ---
        header_match = re.match(r'--- t=([\d\.]+) iter=(\d+) Gen_(\d+) ---', line)
        if header_match:
            time = float(header_match.group(1))
            iteration = int(header_match.group(2))
            gen_id = int(header_match.group(3))
            
            data['times'].append(time)
            data['iterations'].append(iteration)
            
            # Parse residuals line
            i += 1
            if i < len(lines) and lines[i].startswith('RESIDUALS:'):
                residual_line = lines[i].strip()
                residuals = parse_residuals(residual_line)
                for key, val in residuals.items():
                    data['residuals'][key].append(val)
            
            # Parse states line
            i += 1
            if i < len(lines) and lines[i].startswith('STATES:'):
                states = parse_states(lines[i].strip())
                for key, val in states.items():
                    data['states'][key].append(val)
            
            # Parse currents line
            i += 1
            if i < len(lines) and lines[i].startswith('CURRENTS:'):
                currents = parse_currents(lines[i].strip())
                for key, val in currents.items():
                    data['currents'][key].append(val)
            
            # Parse voltages line
            i += 1
            if i < len(lines) and lines[i].startswith('VOLTAGES:'):
                voltages = parse_voltages(lines[i].strip())
                for key, val in voltages.items():
                    data['voltages'][key].append(val)
            
            # Parse power line
            i += 1
            if i < len(lines) and lines[i].startswith('POWER:'):
                power = parse_power(lines[i].strip())
                for key, val in power.items():
                    data['power'][key].append(val)
            
            # Parse equation components
            for _ in range(4):  # Ed'_EQN, Eq'_EQN, SLIP_EQN, AVR_EQN
                i += 1
                if i < len(lines):
                    comp = parse_equation_components(lines[i].strip())
                    for key, val in comp.items():
                        if key in data['equation_components']:
                            data['equation_components'][key].append(val)
            
            # Parse max residual
            i += 1
            if i < len(lines) and lines[i].startswith('MAX_RESIDUAL:'):
                max_res = float(lines[i].split(':')[1].strip())
                data['max_residuals'].append(max_res)
        
        i += 1
    
    return data

def parse_residuals(line):
    """Parse residuals line: RESIDUALS: Ed'=+1.234e-05 Eq'=-2.345e-06 ..."""
    residuals = {}
    parts = line.split()[1:]  # Skip 'RESIDUALS:'
    for part in parts:
        if '=' in part:
            key, val = part.split('=')
            key = key.replace("'", "")  # Remove quotes for dict keys
            residuals[key] = float(val)
    return residuals

def parse_states(line):
    """Parse states line: STATES: Ed'=0.1234 Eq'=0.5678 ..."""
    states = {}
    pattern = r"(\w+'?)=([-+]?[\d\.]+)"
    matches = re.findall(pattern, line)
    for key, val in matches:
        key = key.replace("'", "")  # Remove quotes for dict keys
        states[key] = float(val)
    return states

def parse_currents(line):
    """Parse currents line: CURRENTS: id=0.1234 iq=0.5678"""
    currents = {}
    pattern = r"(\w+)=([-+]?[\d\.]+)"
    matches = re.findall(pattern, line)
    for key, val in matches:
        currents[key] = float(val)
    return currents

def parse_voltages(line):
    """Parse voltages line: VOLTAGES: vd=0.1234 vq=0.5678 ..."""
    voltages = {}
    pattern = r"(\w+)=([-+]?[\d\.]+)"
    matches = re.findall(pattern, line)
    for key, val in matches:
        voltages[key] = float(val)
    return voltages

def parse_power(line):
    """Parse power line: POWER: Pm=0.1234 Pe=0.5678"""
    power = {}
    pattern = r"(\w+)=([-+]?[\d\.]+)"
    matches = re.findall(pattern, line)
    for key, val in matches:
        power[key] = float(val)
    return power

def parse_equation_components(line):
    """Parse equation component lines"""
    components = {}
    pattern = r"(\w+)=([-+]?[\d\.e\-+]+)"
    matches = re.findall(pattern, line)
    for key, val in matches:
        try:
            components[key] = float(val)
        except ValueError:
            pass  # Skip non-numeric values
    return components

def analyze_convergence_issues(data, gen_id):
    """Analyze convergence issues for a specific generator"""
    print(f"\n🔍 CONVERGENCE ANALYSIS FOR GENERATOR {gen_id}")
    print("=" * 60)
    
    if not data['times']:
        print("❌ No data found in log file")
        return
    
    # Find problematic residuals
    max_residuals = np.array(data['max_residuals'])
    problem_indices = np.where(max_residuals > 1e-3)[0]
    
    print(f"📊 Total Newton iterations logged: {len(data['times'])}")
    print(f"⚠️  Iterations with large residuals (>1e-3): {len(problem_indices)}")
    
    if len(problem_indices) > 0:
        print(f"📈 Max residual encountered: {max_residuals.max():.2e}")
        print(f"📉 Min residual achieved: {max_residuals.min():.2e}")
        
        # Find worst residual equation at first problematic iteration
        worst_idx = problem_indices[0]
        worst_residuals = {}
        for eq_name in ['Ed', 'Eq', "Ed''", "Eq''", 'delta', 'slip', 'Efd']:
            if len(data['residuals'][eq_name]) > worst_idx:
                worst_residuals[eq_name] = abs(data['residuals'][eq_name][worst_idx])
        
        worst_equation = max(worst_residuals, key=worst_residuals.get)
        print(f"🎯 Worst equation at first problem: {worst_equation} = {worst_residuals[worst_equation]:.2e}")
        
        # Analyze state evolution
        print(f"\n📋 STATE VALUES AT FIRST PROBLEM (t={data['times'][worst_idx]:.6f}):")
        for state_name in ['Ed', 'Eq', "Ed''", "Eq''", 'delta', 'slip', 'Efd']:
            if len(data['states'][state_name]) > worst_idx:
                val = data['states'][state_name][worst_idx]
                print(f"   {state_name:8s}: {val:+.6f}")
        
        # Analyze equation components
        if len(data['equation_components']['Ed_time']) > worst_idx:
            print(f"\n🔧 EQUATION COMPONENTS AT FIRST PROBLEM:")
            comp = data['equation_components']
            print(f"   Ed' equation  - Time: {comp['Ed_time'][worst_idx]:+.3e}, Steady: {comp['Ed_steady'][worst_idx]:+.3e}")
            print(f"   Eq' equation  - Time: {comp['Eq_time'][worst_idx]:+.3e}, Steady: {comp['Eq_steady'][worst_idx]:+.3e}")
            print(f"   Slip equation - Time: {comp['slip_time'][worst_idx]:+.3e}, Torque: {comp['torque_imb'][worst_idx]:+.3e}")
            if len(comp['Efd_time']) > worst_idx:
                print(f"   AVR equation  - Time: {comp['Efd_time'][worst_idx]:+.3e}, Steady: {comp['Efd_steady'][worst_idx]:+.3e}")
    
    # Check for numerical issues
    print(f"\n⚕️  NUMERICAL HEALTH CHECK:")
    check_numerical_issues(data)

def check_numerical_issues(data):
    """Check for common numerical issues"""
    issues = []
    
    # Check for extremely large state values
    for state_name, values in data['states'].items():
        if values:
            max_val = max(abs(v) for v in values)
            if max_val > 100:
                issues.append(f"Large {state_name} values (max: {max_val:.2e})")
    
    # Check for extremely large currents
    for curr_name, values in data['currents'].items():
        if values:
            max_val = max(abs(v) for v in values)
            if max_val > 10:
                issues.append(f"Large {curr_name} currents (max: {max_val:.2e})")
    
    # Check power balance
    if data['power']['Pm'] and data['power']['Pe']:
        power_imbalances = [abs(pm - pe) for pm, pe in zip(data['power']['Pm'], data['power']['Pe'])]
        max_imbalance = max(power_imbalances) if power_imbalances else 0
        if max_imbalance > 1.0:
            issues.append(f"Large power imbalance (max: {max_imbalance:.2e})")
    
    if issues:
        for issue in issues:
            print(f"   ⚠️  {issue}")
    else:
        print(f"   ✅ No obvious numerical issues detected")

def plot_convergence_history(data, gen_id):
    """Plot convergence history for a generator"""
    if not data['times']:
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Generator {gen_id} Jacobian Solver Debug', fontsize=16)
    
    # Plot residuals
    axes[0,0].semilogy(data['max_residuals'], 'r-o', markersize=3)
    axes[0,0].set_title('Max Residual per Iteration')
    axes[0,0].set_ylabel('Max Residual')
    axes[0,0].grid(True)
    
    # Plot individual residuals
    residual_names = ['Ed', 'Eq', 'slip', 'Efd']
    for name in residual_names:
        if data['residuals'][name]:
            axes[0,1].semilogy(np.abs(data['residuals'][name]), label=name)
    axes[0,1].set_title('Individual Residuals')
    axes[0,1].legend()
    axes[0,1].grid(True)
    
    # Plot state evolution
    state_names = ['delta', 'slip', 'Efd']
    for name in state_names:
        if data['states'][name]:
            axes[0,2].plot(data['states'][name], label=name)
    axes[0,2].set_title('Key State Evolution')
    axes[0,2].legend()
    axes[0,2].grid(True)
    
    # Plot power
    if data['power']['Pm'] and data['power']['Pe']:
        axes[1,0].plot(data['power']['Pm'], 'b-', label='Pm')
        axes[1,0].plot(data['power']['Pe'], 'r-', label='Pe')
        axes[1,0].set_title('Power')
        axes[1,0].legend()
        axes[1,0].grid(True)
    
    # Plot currents
    if data['currents']['id'] and data['currents']['iq']:
        axes[1,1].plot(data['currents']['id'], 'b-', label='id')
        axes[1,1].plot(data['currents']['iq'], 'r-', label='iq')
        axes[1,1].set_title('Currents')
        axes[1,1].legend()
        axes[1,1].grid(True)
    
    # Plot voltages
    voltage_names = ['Vt', 'Vref']
    for name in voltage_names:
        if data['voltages'][name]:
            axes[1,2].plot(data['voltages'][name], label=name)
    axes[1,2].set_title('Terminal Voltage')
    axes[1,2].legend()
    axes[1,2].grid(True)
    
    plt.tight_layout()
    plt.savefig(f'sim/jacobian_debug_gen{gen_id}.png', dpi=150, bbox_inches='tight')
    print(f"📊 Plot saved: sim/jacobian_debug_gen{gen_id}.png")

def main():
    """Main analysis function"""
    print("🔍 JACOBIAN SOLVER DEBUG ANALYZER")
    print("=" * 50)
    
    sim_dir = Path("../sim")
    if not sim_dir.exists():
        print("❌ sim/ directory not found. Run the simulation first.")
        return
    
    # Find all generator log files
    log_files = list(sim_dir.glob("jacobian_gen*.log"))
    
    if not log_files:
        print("❌ No Jacobian generator log files found.")
        print("   Expected files like: sim/jacobian_gen0.log, sim/jacobian_gen1.log, ...")
        return
    
    print(f"📁 Found {len(log_files)} generator log files")
    
    # Analyze each generator
    for log_file in sorted(log_files):
        gen_match = re.search(r'jacobian_gen(\d+)\.log', log_file.name)
        if gen_match:
            gen_id = int(gen_match.group(1))
            print(f"\n🔧 Analyzing {log_file.name}...")
            
            data = parse_generator_log(log_file)
            analyze_convergence_issues(data, gen_id)
            
            # Create plots if matplotlib is available
            try:
                plot_convergence_history(data, gen_id)
            except Exception as e:
                print(f"⚠️  Could not create plots: {e}")
    
    # Check main jacobian log
    main_log = sim_dir / "jacobian.log"
    if main_log.exists():
        print(f"\n📋 Main Jacobian log summary:")
        with open(main_log, 'r') as f:
            lines = f.readlines()
        
        error_lines = [line for line in lines if '❌' in line or 'failed' in line.lower()]
        if error_lines:
            print("🚨 CRITICAL ERRORS FOUND:")
            for line in error_lines[:5]:  # Show first 5 errors
                print(f"   {line.strip()}")
        else:
            print("✅ No critical errors in main log")

if __name__ == "__main__":
    main() 