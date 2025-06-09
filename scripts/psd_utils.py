import os
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Dict, Union

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

# Added helper functions for event timing -------------------------------------
_EVENTS_PATH = Path('sim') / 'events.csv'

def load_events_csv() -> pd.DataFrame:
    """Return dataframe of the events log (time,event,detail). If the file does
    not exist, an empty dataframe is returned to allow graceful fallback.
    """
    if not _EVENTS_PATH.exists():
        # Return empty dataframe with expected columns for downstream safety
        return pd.DataFrame(columns=['time', 'event', 'detail'])
    return pd.read_csv(_EVENTS_PATH)


def get_event_times(event_name: str) -> List[float]:
    """Return a list of times (as float seconds) where *event_name* occurred.
    If the events file is missing or the specific event is not present, an
    empty list is returned. This allows analysis scripts to fall back to their
    heuristic detection methods seamlessly.
    """
    df_evt = load_events_csv()
    if df_evt.empty:
        return []
    times = df_evt.loc[df_evt['event'] == event_name, 'time'].astype(float).tolist()
    return times 