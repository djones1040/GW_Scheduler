from enum import Enum

import numpy as np
from astroplan import Observer
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.coordinates import SkyCoord
from dateutil.parser import parse

class TargetType(Enum):
	Supernova = 1
	Template = 2
	Standard = 3

class Target:
	def __init__(self, name, coord, net_prob, timezone,
				 observatory_lat, observatory_lon, observatory_elev, sidereal_radian_array, \
				 disc_date=None, apparent_mag=None, obs_date=None, obs_datetime=None,
				 obs_filter=None, gw_exptime=None, tstart=None):
		# Provided by Constructor
		self.name = name
		self.coord = coord
		self.net_prob = net_prob
		self.disc_date = disc_date
		self.apparent_mag = apparent_mag
		self.obs_date = obs_date
		self.obs_datetime = obs_datetime
		self.observatory_lat = observatory_lat
		self.observatory_lon = observatory_lon
		self.observatory_elev = observatory_elev
		self.timezone = timezone
		self.gw_exptime = gw_exptime
		self.obs_filter = obs_filter
		self.tstart = tstart
		
		# Computed by Constructor
		#self.raw_airmass_array = self.compute_airmass(observatory_lat, sidereal_radian_array)
		
		# Computed by Telescope
		self.starting_index = 0 # Used to order net priority
		self.exposures = None # Dictionary: filter:minutes
		self.total_observable_min = 0 # How many minutes in the night is the target observable?
		self.total_minutes = 0 # Total length of observation
		self.fraction_time_obs = 9999 # TotalMinutes / TotalObservableMin
		self.total_good_air_mass = 9999 # Proxy for elevation
		self.scheduled_time_array = None # Airmass plot abscissa
		self.scheduled_airmass_array = None # Airmass plot ordinate
	
	def compute_airmass(self, observatory_lon, observatory_lat, observatory_height, timezone,
						current_time):

		location = EarthLocation.from_geodetic(observatory_lon*u.deg, observatory_lat*u.deg,
											   observatory_height*u.m)
		telescope = Observer(location=location, timezone=timezone)

		target_coord = SkyCoord(self.coord.ra.deg,self.coord.dec.deg,unit=u.deg)
		airmass = telescope.altaz(current_time, target_coord).secz
		
		return airmass

	def compute_rise_time(self, observatory_lon, observatory_lat, observatory_height, timezone,
						  current_time, airmass_at_rise=3):

		location = EarthLocation.from_geodetic(observatory_lon*u.deg, observatory_lat*u.deg,
											   observatory_height*u.m)
		telescope = Observer(location=location, timezone=timezone)

		target_coord = SkyCoord(self.coord.ra.deg,self.coord.dec.deg,unit=u.deg)

		horizon = 90-np.arccos(1/airmass_at_rise)*(180./np.pi)
		rise_time = telescope.target_rise_time(current_time, target_coord, horizon=horizon*u.deg, which='next')
		
		return parse(rise_time.isot)

	def compute_set_time(self, observatory_lon, observatory_lat, observatory_height, timezone,
						 current_time, airmass_at_set=3):

		location = EarthLocation.from_geodetic(observatory_lon*u.deg, observatory_lat*u.deg,
											   observatory_height*u.m)
		telescope = Observer(location=location, timezone=timezone)

		target_coord = SkyCoord(self.coord.ra.deg,self.coord.dec.deg,unit=u.deg)

		horizon = 90-np.arccos(1/airmass_at_set)*(180./np.pi)
		set_time = telescope.target_set_time(current_time, target_coord, horizon=horizon*u.deg, which='next')
		
		return parse(set_time.isot)
