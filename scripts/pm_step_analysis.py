import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from psd_utils import load_generator_csv, save_fig, nice_axes, get_event_times

def main():
    df = load_generator_csv(0)
    time = df['time']
    pm = df['mech_power']
    # Prefer explicit event log when available
    evt_times = get_event_times('pm_step')
    if evt_times:
        step_time = evt_times[0]
        # find closest index for annotation
        step_idx = (df['time'] - step_time).abs().idxmin()
    else:
        dpm = pm.diff().abs()
        step_idx = dpm.idxmax() if dpm.max() > 1e-4 else 0
        step_time = time.iloc[step_idx]

    fig, ax = plt.subplots(3,1, figsize=(8,8), sharex=True)
    ax[0].plot(time, pm, label='Pm', color='tab:blue'); ax[0].axvline(step_time, ls='--'); ax[0].set_ylabel('Pm (pu)'); nice_axes(ax[0])
    ax[1].plot(time, df['elec_power'], label='Pe', color='tab:red'); ax[1].axvline(step_time, ls='--'); ax[1].set_ylabel('Pe (pu)'); nice_axes(ax[1])
    ax[2].plot(time, df['slip'], label='Slip', color='tab:green'); ax[2].axvline(step_time, ls='--'); ax[2].set_ylabel('Slip (pu)'); nice_axes(ax[2]); ax[2].set_xlabel('Time (s)')

    save_fig(fig, 'pm_time_series.png')

    summary = Path('report_plot') / 'pm_summary.txt'
    with summary.open('w') as fp:
        fp.write(f"Pm step at t = {step_time:.2f} s\n")
        fp.write(f"ΔPm = {pm.iloc[step_idx] - pm.iloc[step_idx-1]:.3f} pu\n")
        fp.write(f"Final Slip = {df['slip'].iloc[-1]:.5f} pu\n")
    print(f"[summary] Saved {summary}")

if __name__=='__main__':
    main() 