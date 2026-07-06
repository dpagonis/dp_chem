import re
from datetime import datetime, timedelta, timezone
import pandas as pd
import warnings

#####################################

def parse_filename(filename):
    match = re.match(r'([^_]+)_([^_]+)_(\d{8})_R([A-Z]|\d+)(_L(\d+))?(.*)?\.ict', filename)
    if match:
        measurement_id = match.group(1)
        platform = match.group(2)
        date = match.group(3)
        revision = match.group(4)
        leg = match.group(6) if match.group(6) else None
        return measurement_id, platform, date, revision, leg
    return None, None, None, None, None

##################################

def sort_revision_key(revision):
    if revision.isdigit():
        return (1, int(revision))
    else:
        return (0, revision)

# Function to get the latest file based on revision
def get_latest_file(files):
    # Sort files by revision
    files.sort(key=lambda x: sort_revision_key(x[0]), reverse=True)
    # Return the latest file (last in sorted list)
    return files[0][1]

####################################

def dt_sec_since_midnight(dt : datetime):
    """return the number of seconds since midnight, UTC
    for naive datetime objects, assumed timezone is UTC"""

    if dt.tzinfo is None:
        dt_utc = dt.replace(tzinfo=timezone.utc)
    else:
        dt_utc = dt.astimezone(timezone.utc)
    s = dt_utc.hour * 3600 + dt_utc.minute * 60 + dt_utc.second + dt_utc.microsecond * 1e-6
    return s

####################################
def load(filepath):
    with open(filepath, 'r') as file:
        # fetch header 
        header_length = int(file.readline().split(',')[0])-1
        header = [file.readline().strip() for _ in range(header_length - 1)]
        
        # Extract the start date from the header (line 7, index 6)
        start_date_line = header[5].split(',')
        start_date = datetime(year=int(start_date_line[0]), month=int(start_date_line[1]), day=int(start_date_line[2]))
        
        # Read the data into a DataFrame (note skiprows is zero since we're doing open & readline)
        data = pd.read_csv(file, skiprows=0)
        
        utc = timezone.utc
        
        # Convert the icartt time column to UTC datetime objects
        # Some people don't follow the icartt standard and instead give iso strings of the UTC timestamp. we handle that, but whine about it

        utc = timezone.utc
        try:
            data['Time_Start'] = data.iloc[:, 0].apply(lambda x: utc.localize(start_date + timedelta(seconds=float(x))))
        except:
            try:
                data['Time_Start'] = data.iloc[:, 0].apply(lambda x: datetime.fromisoformat(x).replace(tzinfo=utc))
                msg = f'file {filepath} does not comply with icartt standard time format. parsing timestamps as ISO-format with assumed timzeone UTC'
                warnings.warn(msg, UserWarning) 
            except:
                raise ValueError(f'file {filepath} does not comply with icartt standard time format and could not be loaded.') 

        # drop extra 'time' columns and set the datetime column as the index
        data = data.loc[:, ~data.columns.str.contains('time', case=False) | data.columns.isin(['Time_Start'])]
        data.set_index('Time_Start', inplace=True)
        
        return data
