import pandas as pd
import matplotlib.pyplot as plt
from psd_utils import load_generator_csv, nice_axes, save_fig
from pathlib import Path


def main():
    df = load_generator_csv(0)
    time = df['time']

    # detect Vref step (first non-zero diff)
    vref = df['Vref']
    dv = vref.diff().abs()
    step_idx = dv.idxmax() if dv.max() > 1e-4 else 0
    step_time = time.iloc[step_idx]

    fig, ax = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

    ax[0].plot(time, vref, label='Vref', color='tab:blue')
    ax[0].axvline(step_time, ls='--', color='k')
    ax[0].set_ylabel('Vref (pu)'); nice_axes(ax[0])

    ax[1].plot(time, df['Efd'], label='Efd', color='tab:green')
    ax[1].axvline(step_time, ls='--', color='k'); ax[1].set_ylabel('Efd (pu)'); nice_axes(ax[1])

    ax[2].plot(time, df['Vt'], label='Vt', color='tab:orange')
    ax[2].axvline(step_time, ls='--', color='k'); ax[2].set_ylabel('Vt (pu)'); nice_axes(ax[2])

    ax[3].plot(time, df['slip'], label='Slip', color='tab:red')
    ax[3].axvline(step_time, ls='--', color='k'); ax[3].set_ylabel('Slip (pu)'); nice_axes(ax[3])
    ax[3].set_xlabel('Time (s)')

    save_fig(fig, 'vref_time_series.png')

    # simple summary
    summary = Path('report_plot') / 'vref_summary.txt'
    with summary.open('w') as fp:
        fp.write(f"Vref step detected at t = {step_time:.2f} s\n")
        fp.write(f"ΔVref = {vref.iloc[step_idx] - vref.iloc[step_idx-1]:.3f} pu\n")
        fp.write(f"Final Efd = {df['Efd'].iloc[-1]:.3f} pu\n")
        fp.write(f"Final Slip = {df['slip'].iloc[-1]:.5f} pu\n")
    print(f"[summary] Saved {summary}")


if __name__ == '__main__':
    main() 