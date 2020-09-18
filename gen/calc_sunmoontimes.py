import astroplan
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astroplan import Observer
from datetime import datetime

def calc_sunmoontimes(Earth_Location,date,timezone='UTC',printit=False):
    # E.g. Earth_Location = "Cerro Tololo"
    # E.g. date = '2015-12-18'
    location   = EarthLocation.of_site(Earth_Location)
    scope      = Observer(location = location, timezone = timezone,
                                 name = "name", description = "scope")
    t_midnight = Time(date)
    sunset     = scope.sun_set_time(t_midnight, which='nearest')
    eve_twil   = scope.twilight_evening_astronomical(t_midnight, which='nearest')
    midnight   = scope.midnight(t_midnight, which='next')
    morn_twil  = scope.twilight_morning_astronomical(t_midnight, which='next')
    sunrise    = scope.sun_rise_time(t_midnight, which='next')
    moonrise   = scope.moon_rise_time(t_midnight, which='next')
    moonset    = scope.moon_set_time(t_midnight, which='next')

    if printit:
        print("Sunset                        {0.iso} UTC".format(sunset))
        print("Astronomical evening twilight {0.iso} UTC".format(eve_twil))
        print("Midnight                      {0.iso} UTC".format(midnight))
        print("Astronomical morning twilight {0.iso} UTC".format(morn_twil))
        print("Sunrise                       {0.iso} UTC".format(sunrise))
        print("Moonrise                      {0.iso} UTC".format(moonrise))
        print("Moonset                       {0.iso} UTC".format(moonset))

    sunrise_datetime   = datetime.strptime('{0.iso}'.format(sunset), "%Y-%m-%d %H:%M:%S.%f")
    sunset_datetime    = datetime.strptime('{0.iso}'.format(sunrise), "%Y-%m-%d %H:%M:%S.%f")

    eve_twil_datetime  = datetime.strptime('{0.iso}'.format(eve_twil), "%Y-%m-%d %H:%M:%S.%f")
    morn_twil_datetime = datetime.strptime('{0.iso}'.format(morn_twil), "%Y-%m-%d %H:%M:%S.%f")

    moonrise_datetime  = datetime.strptime('{0.iso}'.format(moonrise), "%Y-%m-%d %H:%M:%S.%f")
    moonset_datetime   = datetime.strptime('{0.iso}'.format(moonset), "%Y-%m-%d %H:%M:%S.%f")

    #moonrise_datetime  = datetime.strptime('2015-12-18 07:30:19.330699', "%Y-%m-%d %H:%M:%S.%f")

    sunmoondict={'moonrise':moonrise_datetime,
                 'moonset':moonset_datetime,
                 'sunrise':sunrise_datetime,
                 'sunset':sunset_datetime,
                 'eve_twil':eve_twil_datetime,
                 'morn_twil':morn_twil_datetime
                }
    
    return sunmoondict
