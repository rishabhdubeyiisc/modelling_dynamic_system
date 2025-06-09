import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

REPORT_DIR = Path('report_plot')
REPORT_DIR.mkdir(exist_ok=True)

def load_generator_csv(gen_idx: int = 0) -> pd.DataFrame:
    """Return dataframe for a given generator CSV located in sim/"""
    csv_path = Path('sim') / f'gen{gen_idx}.csv'
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV {csv_path} not found. Run simulation first.")
    df = pd.read_csv(csv_path)
    return df

def nice_axes(ax):
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Time (s)')

def save_fig(fig, name: str):
    fig.tight_layout()
    out = REPORT_DIR / name
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[plot] Saved {out}") 