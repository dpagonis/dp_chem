import inspect
import os
import numpy as np
import getpass
from datetime import datetime, timedelta, timezone
import math
from zoneinfo import ZoneInfo
import math
import numpy as np

from dp_chem import molecule
from dp_chem.units import units
from dp_chem import timeseries


EARTHRADIUS_m = 6.371e6
LOSCHMIDT_mQcm3 = 2.6867811e19
LOSCHMIDT_mQm3 = 2.6867811e25
AVOGADRO = 6.02214076e23


def convert(value:float, unit:str, convert_to:str):
    base_unit = units(unit)
    conversion_factor = base_unit.convert_to(convert_to)
    
    if conversion_factor is None:
        raise ValueError(f"dp_chem.units could not convert {unit} to {convert_to}")
    
    return value * conversion_factor

def dt64_to_dt(dt64:np.datetime64, tz=None):
    dt = dt64.astype('M8[us]').item()
    dt_utc = dt.replace(tzinfo=timezone.utc)
    
    if tz is None:
        return dt_utc
    else:    
        if isinstance(tz,(timezone,ZoneInfo)):
            dt_tz = dt_utc.astimezone(tz)
        elif isinstance(tz,str):
            dt_tz = dt_utc.astimezone(ZoneInfo(tz))
        else:
            raise TypeError(f"dt64_to_dt cannot accept {type(tz)} for tz. Use ZoneInfo or a string")
        return dt_tz

        
def save_fig(fig, fname=None, dpi=200, annotate=True, fontsize=8):
    """
    Save a matplotlib.figure.Figure to disk. If fname is None, a default path
    is produced from main_dir() with a .png extension. Optionally annotate the
    figure with the user, save datetime, and the script that created it.
    """
    import matplotlib.pyplot as _plt

    # determine filename
    if fname is None:
        fname = fig_dir()

    # annotation
    if annotate:
        try:
            user = getpass.getuser()
        except Exception:
            user = os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"

        try:
            main_path = main_dir()
            script_name = os.path.basename(main_path)
        except Exception:
            script_name = "unknown"

        date_str = datetime.now().strftime("%Y-%m-%d %H:%M")
        annotation = f"{user} | {date_str} | {script_name}"

        # Add annotation in the lower-right corner of the figure
        # Use fig.text so it appears on all subplots and is saved with the figure.
        fig.text(0.99, 0, annotation, ha="right", va="bottom",
                 fontsize=fontsize, alpha=0.5)

    # save
    # Accept either a figure object or the pyplot state
    if hasattr(fig, "savefig"):
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    else:
        # assume pyplot
        _plt.savefig(fname, dpi=dpi, bbox_inches="tight")

    return fname

def fig_dir():
    maindir = main_dir()
    return maindir.replace(".py", ".png")

def main_dir():
    # Find the outermost script (__main__)
    for frame in inspect.stack():
        module = inspect.getmodule(frame.frame)
        if module and module.__name__ == "__main__":
            parent_script = getattr(module, "__file__", None)
            if parent_script:
                break
    else:
        raise RuntimeError("Could not determine the outermost script (__main__).")

    parent_dir = os.path.dirname(os.path.abspath(parent_script))
    script_name = os.path.basename(parent_script)
    main_path = os.path.join(parent_dir, script_name)
    return main_path

def MW(s):
    if isinstance(s,str):
        m=molecule(s)
    elif isinstance(s,molecule):
        m=s
    else:
        raise TypeError(f"pass MW a molecule or a string. {type(s)}")
    return f'{m.formula}: {m.molecular_weight.as_num()} g/mol'

def N(T_K,P_mbar):
    mol_m3 = P_mbar*100/(8.3145*T_K)
    molec_cm3 = mol_m3 * 1e-6 * 6.02214076e23
    return molec_cm3

def distance(ll_tuple1,ll_tuple2):
    lat1, lon1 = ll_tuple1
    lat2, lon2 = ll_tuple2
    lat1_r= lat1 *math.pi /180
    lon1_r= lon1 *math.pi /180
    lat2_r= lat2 *math.pi /180
    lon2_r= lon2 *math.pi /180

    d_r = 2*math.asin(math.sqrt((math.sin((lat1_r-lat2_r)/2))**2 + math.cos(lat1_r)*math.cos(lat2_r)*(math.sin((lon1_r-lon2_r)/2))**2))
    d_m = d_r * EARTHRADIUS_m
    return d_m

def mda8(ts):
    """
    calculates maximum daily 8-hour average for every day in a timeseries object
    returns a timeseries object with mda8 values
    """

    # Storage for unique days and their max 8-hour
    unique_days = []
    max_8h_values = []

    # Loop through unique days
    for day in np.unique([t.date() for t in ts.t]):  
        daily_mask = np.array([t.date() == day for t in ts.t])  # Mask for the day's data
        daily_times = ts.t[daily_mask]
        daily_values = ts.data[daily_mask]
        
        max_8h_avg = None
        
        
        # Slide an 8-hour window across the day's data
        for i, center_time in enumerate(daily_times):
            start_time = center_time + timedelta(hours=-4)
            end_time = center_time + timedelta(hours=4)
            mask = (daily_times >= start_time) & (daily_times < end_time)
            
            if np.sum(mask) > 0:  # Ensure there's data in the window
                avg_8h = np.nanmean(daily_values[mask])
                if max_8h_avg is None or avg_8h > max_8h_avg:
                    max_8h_avg = avg_8h
                    dt_day = center_time.date()
                    dt = datetime(dt_day.year,dt_day.month,dt_day.day,0,0,0,tzinfo=center_time.tzinfo)

        # Store results
        unique_days.append(dt)
        max_8h_values.append(max_8h_avg)


    unique_days = np.array(unique_days)
    max_8h_values = np.array(max_8h_values)

    MDA8=timeseries(unique_days,max_8h_values)
    return MDA8

def sza(lat, lon, time):
    from dp_chem import solar
    zenithangle, _ = solar.solar_calcs(lat,lon,time)
    return zenithangle

if __name__ == '__main__':
    print(convert(4253, 'ft','km'))