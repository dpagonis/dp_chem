import pandas as pd 
from pathlib import Path
from zoneinfo import ZoneInfo

# load in the locations dataset
_this_file_dir = Path(__file__).resolve().parent
_csvpath = Path.joinpath(_this_file_dir,'tables/','locations.csv')
_DF = pd.read_csv(_csvpath)

# build the coordinates dictionary
def _build_coords_dict(df):
    latlon_dict = {}
    for i, row in df.iterrows():
        name = row['name']
        aka = row['aka']
        lat	= row['lat']
        lon = row['lon']

        latlon_dict[name] = (lat,lon)

        if isinstance(aka,str):
            aka_list = aka.split(',')
            for alias in aka_list:
                key = alias.strip()
                latlon_dict[key] = (lat,lon)

    return latlon_dict

COORDINATES = _build_coords_dict(_DF)

##########################

def _get_row(site) -> pd.Series:
    for i, row in _DF.iterrows():
        if isinstance(row['aka'],str):
            name_list = [s.strip() for s in row['aka'].split(',')]
            name_list.append(row['name'])
        else:
            name_list = [row['name']]
        
        if site in name_list:
            return row
        
    raise ValueError(f"dp_chem.locations did not find '{site}'")

##########################

def altitude_km(site:str) -> float:
    row = _get_row(site)
    altitude = float(row['altitude_km'])
    return altitude

def altitude_m(site:str) -> float:
    return 1000 * altitude_km(site)

#########################

def timezone(site:str) -> ZoneInfo:
    row = _get_row(site)
    tz_str = row['tz']
    return ZoneInfo(tz_str)

#########################

if __name__ == "__main__":
    print(timezone('ogden'))