from enum import Enum

import numpy as np

class TargetType(Enum):
    Supernova = 1
    Template = 2
    Standard = 3

class Target:
    def __init__(self, name, coord, priority, target_type, observatory_lat, sidereal_radian_array, \
                 disc_date=None, apparent_mag=None, obs_date=None):
        # Provided by Constructor
        self.name = name
        self.coord = coord
        self.priority = priority
        self.type = target_type
        self.disc_date = disc_date
        self.apparent_mag = apparent_mag
        self.obs_date = obs_date
        
        # Computed by Constructor
        self.raw_airmass_array = self.compute_airmass(observatory_lat, sidereal_radian_array)
        
        # Computed by Telescope
        self.net_priority = self.priority
        self.starting_index = 0 # Used to order net priority
        self.exposures = None # Dictionary: filter:minutes
        self.total_observable_min = 0 # How many minutes in the night is the target observable?
        self.total_minutes = 0 # Total length of observation
        self.fraction_time_obs = 9999 # TotalMinutes / TotalObservableMin
        self.total_good_air_mass = 9999 # Proxy for elevation
        self.scheduled_time_array = None # Airmass plot abscissa
        self.scheduled_airmass_array = None # Airmass plot ordinate
    
    def compute_airmass(self, observatory_lat, sidereal_radian_array):
        n = len(sidereal_radian_array)

        RA = np.empty(n)
        RA.fill(self.coord.ra.radian)

        DEC = np.empty(n)
        DEC.fill(self.coord.dec.radian)

        HA = sidereal_radian_array - RA
        LAT = np.empty(n)
        LAT.fill(observatory_lat)
        
        term1 = np.sin(DEC)*np.sin(LAT)
        term2 = np.cos(DEC)*np.cos(LAT)*np.cos(HA)
        am = 1.0/(np.sin(np.arcsin(term1+term2)))
        
        am[(am > 3.0) | (am < 1.0)] = 9999

        return am