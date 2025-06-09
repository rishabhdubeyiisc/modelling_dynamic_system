import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from psd_utils import load_generator_csv, save_fig, nice_axes, get_event_times

def main():
    df = load_generator_csv(0)
    time = df['time']
    vq = df['VQ'] if 'VQ' in df else df['vq']
    # Prefer explicit fault timing from events.csv
    evt_start = get_event_times('fault_start')
    evt_end   = get_event_times('fault_end')

    if evt_start and evt_end:
        start_time = evt_start[0]
        clear_time = evt_end[0]
        start_idx = (df['time'] - start_time).abs().idxmin()
        clear_idx = (df['time'] - clear_time).abs().idxmin()
    else:
        # fallback heuristic: big drop in VQ indicates fault start
        dvq = vq.diff()
        start_idx = dvq.idxmin()
        start_time = time.iloc[start_idx]
        # fault clear when dvq largest positive after start
        clear_idx = dvq[start_idx+1:].idxmax()
        clear_time = time.iloc[clear_idx]

    fig, ax = plt.subplots(2,1, figsize=(8,6), sharex=True)
    ax[0].plot(time, vq, label='VQ'); ax[0].axvline(start_time, ls='--', color='r'); ax[0].axvline(clear_time, ls='--', color='g'); ax[0].set_ylabel('VQ (pu)'); nice_axes(ax[0])
    ax[1].plot(time, df['delta'], label='delta'); ax[1].axvline(start_time, ls='--', color='r'); ax[1].axvline(clear_time, ls='--', color='g'); ax[1].set_ylabel('δ (rad)'); ax[1].set_xlabel('Time (s)'); nice_axes(ax[1])

    save_fig(fig, 'fault_time_series.png')

    summary = Path('report_plot') / 'fault_stability_report.txt'
    delta_max = df['delta'].max(); delta_min = df['delta'].min()
    with summary.open('w') as fp:
        fp.write(f"Fault start at t = {start_time:.2f} s, clear at t = {clear_time:.2f} s\n")
        fp.write(f"Max δ = {delta_max:.2f} rad, Min δ = {delta_min:.2f} rad\n")
        fp.write("Rotor angle stable" if abs(delta_max - delta_min) < 3.14 else "Rotor angle unstable")
    print(f"[summary] Saved {summary}")

if __name__=='__main__':
    main() 