from datetime import timezone
import math

def sza(lat, lon, time):
    zenithangle, _ = solar_calcs(lat,lon,time)
    return zenithangle

#############################

def _declin_eot(time):

    if time.tzinfo is None: # UTC if we are tz-aware
        time = time.replace(tzinfo=timezone.utc)

    julian_day = time.timestamp() / 86400.0 + 2440587.5
    julian_century =(julian_day-2451545)/36525
    geom_mean_lon_sun = (280.46646 + julian_century * (36000.76983 + julian_century * 0.0003032)) % 360
    geom_mean_anomaly_sun =357.52911+julian_century*(35999.05029 - 0.0001537*julian_century)
    geom_mean_anomaly_sun_rad = geom_mean_anomaly_sun*math.pi/180
    eccent_earth_orbit = 0.016708634-julian_century*(0.000042037+0.0000001267*julian_century)
    sun_eqn_center =math.sin(geom_mean_anomaly_sun_rad)*(1.914602-julian_century*(0.004817+0.000014*julian_century))+math.sin(2*geom_mean_anomaly_sun_rad)*(0.019993-0.000101*julian_century)+math.sin(3*geom_mean_anomaly_sun_rad)*0.000289
    sun_true_long = geom_mean_lon_sun + sun_eqn_center
    sun_true_anom = geom_mean_anomaly_sun + sun_eqn_center
    sun_rad_vector_AUs =(1.000001018*(1-eccent_earth_orbit*eccent_earth_orbit))/(1+eccent_earth_orbit*math.cos(sun_true_anom*math.pi/180))   
    sun_app_long = sun_true_long-0.00569-0.00478*math.sin((125.04-1934.136*julian_century)*math.pi/180)
    mean_obliq_ecliptic = 23+(26+((21.448-julian_century*(46.815+julian_century*(0.00059-julian_century*0.001813))))/60)/60
    obliq_corr = mean_obliq_ecliptic+0.00256*math.cos((125.04-1934.136*julian_century)*math.pi/180)
    sun_rt_ascen = math.degrees(math.atan2(math.cos(math.radians(obliq_corr)) * math.sin(math.radians(sun_app_long)), math.cos(math.radians(sun_app_long))))
    sun_declin = math.degrees(math.asin(math.sin(math.radians(obliq_corr))*math.sin(math.radians(sun_app_long))))
    var_y = math.tan(math.radians(obliq_corr/2))*math.tan(math.radians(obliq_corr/2))
    eq_of_time =4*math.degrees(var_y*math.sin(2*math.radians(geom_mean_lon_sun))-2*eccent_earth_orbit*math.sin(math.radians(geom_mean_anomaly_sun))+4*eccent_earth_orbit*var_y*math.sin(math.radians(geom_mean_anomaly_sun))*math.cos(2*math.radians(geom_mean_lon_sun))-0.5*var_y*var_y*math.sin(4*math.radians(geom_mean_lon_sun))-1.25*eccent_earth_orbit*eccent_earth_orbit*math.sin(2*math.radians(geom_mean_anomaly_sun)))

    return sun_declin, eq_of_time

#############################

def solar_calcs(lat,lon,time):
    """
    Calculate:
    1. solar zenith angle in degrees
    2. solar azimuth, degrees cw from N  
    
    Coords northing, easting. datetime object for time. Be tz aware, or UTC will be assumed

    The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators are based on equations from Astronomical Algorithms, by Jean Meeus. 
    The sunrise and sunset results are theoretically accurate to within a minute for locations between +/- 72° latitude, and within 10 minutes 
    outside of those latitudes. However, due to variations in atmospheric composition, temperature, pressure and conditions, observed values may 
    vary from calculations.
    
    Rewritten from https://gml.noaa.gov/grad/solcalc/calcdetails.html
    """
    if time.tzinfo is None: # assume UTC if we are tz-naive
        time = time.replace(tzinfo=timezone.utc)
    
    sun_declin, eq_of_time = _declin_eot(time)
   
    utc_offset = time.utcoffset().total_seconds() / 3600  # Get UTC offset in hours
    time_past_midnight = time.hour * 60 + time.minute + time.second / 60  # Time in minutes past midnight
    true_solar_time = (time_past_midnight + eq_of_time + 4 * lon - 60 * utc_offset) % 1440
    hour_angle = true_solar_time/4+180 if true_solar_time<0 else true_solar_time/4-180
    solar_zenith_angle =math.degrees(math.acos(math.sin(math.radians(lat))*math.sin(math.radians(sun_declin))+math.cos(math.radians(lat))*math.cos(math.radians(sun_declin))*math.cos(math.radians(hour_angle))))
    if hour_angle > 0:
        azimuth = (math.degrees(math.acos(((math.sin(math.radians(lat)) * math.cos(math.radians(solar_zenith_angle))) - math.sin(math.radians(sun_declin))) / (math.cos(math.radians(lat)) * math.sin(math.radians(solar_zenith_angle))))) + 180) % 360
    else:
        azimuth = (540 - math.degrees(math.acos(((math.sin(math.radians(lat)) * math.cos(math.radians(solar_zenith_angle))) - math.sin(math.radians(sun_declin))) / (math.cos(math.radians(lat)) * math.sin(math.radians(solar_zenith_angle)))))) % 360
    return solar_zenith_angle, azimuth