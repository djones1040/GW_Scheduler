import configparser
import os

import Constants
from Observatory import Observatory
from Telescope import Swope, Nickel
from Utilities import *
from Target import TargetType, Target

from dateutil.parser import parse
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as unit

from datetime import datetime, timedelta

class Scheduler:
	def __init__(self):
		self.telescopes = []

	def addwarning(self,warning):
		print(warning)
		self.warnings.append(warning)
		
	def add_options(self, parser=None, usage=None, config=None):
		if parser == None:
			parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")
			
		parser.add_argument("-f", "--file", default=config.get('main','file'),
							help="CSV file with targets to schedule.")
		parser.add_argument("-d", "--date", default=config.get('main','file'),
							help="YYYYMMDD formatted observation date.")
		parser.add_argument("--total_prob_goal", default=config.get('main','total_prob_goal'),
							help="total probability goal.  Helps to avoid low-prob targets with my simplified weighting scheme")
		parser.add_argument("--airmass_max", default=config.get('main','airmass_max'),
							help="maximum airmass to observe at")
		
		for telescope in config.sections():
			if telescope == 'main': continue
			self.telescopes += [telescope]
			
			parser.add_argument("--%s_slewtime_per_deg",default=config.get(telescope,'slewtime_per_deg'))
			parser.add_argument("--%s_plate_scale",default=config.get(telescope,'plate_scale'))
			parser.add_argument("--%s_detector_width_npix",default=config.get(telescope,'detector_width_npix'))
			parser.add_argument("--%s_detector_height_npix",default=config.get(telescope,'detector_height_npix'))
			parser.add_argument("--%s_lat",default=config.get(telescope,'lat'))
			parser.add_argument("--%s_lon",default=config.get(telescope,'lon'))
			parser.add_argument("--%s_elevation",default=config.get(telescope,'elevation'))
			parser.add_argument("--%s_utc_offset",default=config.get(telescope,'utc_offset'))
			parser.add_argument("--%s_use",default=config.get(telescope,'use'),help="schedule observations on this telescope?")
			
		return parser

	def main(self):
	
		file_name = self.args.file
		obs_date = self.args.date
		preview_plot = self.args.plot

		observatories = {}
		targets_to_schedule = {}
		fov = []
		for t in self.telescopes:
			tel = Observatory(
				name=t,
				lon=self.args.__dict__['%s_lon'%t],
				lat=self.args.__dict__['%s_lat'%t],
				elevation=self.args.__dict__['%s_elevation'%t],
				horizon="-12",
				telescopes={"Swope":Swope()},
				obs_date_str=obs_date,
				utc_offset=self.args.__dict__['%s_utc_offset'%t],
				fov=self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]/60.,
				gw_exptime=self.args.__dict__['%s_gw_exptime'%t],
				slewtime_per_deg=self.args.__dict__['%s_slew_time_per_deg'%t]
			)
			observatories[t] = tel
			targets_to_schedule[t] = []
			fov += [self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]/60.]
			
		target_data = get_targets("%s" % file_name)
		names = [t[0] for t in target_data]
		ra = [t[1] for t in target_data]
		dec = [t[2] for t in target_data]
		priorities = [float(t[3]) for t in target_data]
		disc_dates = [t[4] for t in target_data]
		disc_mags = [float(t[5]) for t in target_data]
		types = [t[6] for t in target_data]
		coords = SkyCoord(ra,dec,unit=(unit.hour, unit.deg))

		# sort by priorities
		iSort = np.argsort(priorities)[::-1]
		names,ra,dec,priorities,disc_dates,disc_mags,types = \
			names[iSort],ra[iSort],dec[iSort],priorities[iSort],\
			disc_dates[iSort],disc_mags[iSort],types[iSort]


		# datetime array
		def datetime_range(start, end, delta):
			current = start
			while current < end:
				yield current
				current += delta

		dts = [dt.strftime(obs_date) for dt in 
			   datetime_range(datetime(int(obsdate[0:4]), int(obsdate[4:6]), int(obsdate[6:8]), 12),
							  datetime(int(obsdate[0:4]), int(obsdate[4:6]), int(obsdate[6:8])+1, 12), 
							  timedelta(minutes=1))]
		telra = [None]*len(self.telescopes)
		teldec = [None]*len(self.telescopes)
		current_obs_tstart = [None]*len(self.telescopes)
		for d in dts:
			
			for obs,telra,teldec,i in zip(
					self.telescopes[np.argsort(fov)],telra[np.argsort(fov)],
					teldec[np.argsort(fov)],range(len(teldec))):

				# keep going if it's not night yet
				if obs.utc_begin_night > d: continue
				if obs.utc_end_night < d: continue

				if d - current_obs_tstart < obs.gw_exptime:
					continue
				
				targetlist,airmass_list,time_to_set_list,slewtime_list = np.array([]),np.array([]),np.array([]),np.array([])

				probsum = 0
				for i in range(len(prob)):
					probsum += prob[i]
					if probsum > self.args.total_prob_goal: iProbMax = i
					
				for n,c,p in zip(names[0:iProbMax],coords[0:iProbMax],prob[0:iProbMax]):
				
					target = Target(
						name=names[j], 
						coord=coords[j], 
						priority=priorities[j], 
						target_type=target_type, 
						observatory_lat=obs.ephemeris.lat, 
						sidereal_radian_array=obs.sidereal_radian_array, 
						disc_date=disc_date, 
						apparent_mag=disc_mags[j], 
						obs_date=obs.obs_date,
						obs_datetime=d
					)
					airmass = target.compute_airmass(target.observatory_lat,obs.sidereal_radian_array)
					time_to_set = d - target.compute_time_until_rise(target.observatory_lon,target.observatory_lat,target.observatory_height,
																	 target.timezone,d,airmass_at_set=self.args.airmass_max)
					time_to_rise = d - target.compute_time_until_set(target.observatory_lon,target.observatory_lat,target.observatory_height,
																	 target.timezone,d,airmass_at_set=self.args.airmass_max)

					if telra and teldec:
						telsc = SkyCoord(telra,teldec,unit=u.deg)
						objsc = SkyCoord(c[0],c[i],unit=u.deg)
						sep_deg = telsc.separation_to(objsc).deg
						slewtime = obs.slewtime_per_deg*sep_deg
					else:
						slewtime = 0
						
					targetlist = np.append(targetlist,target)
					airmass_list = np.append(airmass_list,airmass)
					time_to_set_list = np.append(time_to_set_list,time_to_set)
					time_to_rise_list = np.append(time_to_rise_list,time_to_rise)
					slewtime_list = np.append(slewtime_list,slewtime)

				# compute modified priorities
				# if rising and AM > 2, lower priority (-1)
				# if setting w/i 30 min (AM > 2.5?), raise priority (+=3)
				# if setting w/i 60 min (AM > 2.5?), raise priority (+=2)
				# if slewtime < 15s, raise priority (+=1)
				# only look at tiles that comprise 50% of the total probability space
				# file should have a "done" flag, so you can get 50% of the total
				#      probability w/o scheduling finished things.  Then you can
				#      increase the threshold to 60%, 70%, etc as the night goes on
				priority = np.ones(len(targetlist))
				priority[(airmass_list > 2) & (time_to_rise_list < 1)] -= 1
				priority[(time_to_set_list < 60)] += 2
				priority[(time_to_set_list < 30)] += 1
				priority[(slewtime_list < 15)] += 1

				# then have to schedule things for x number of minutes
				if i == 0:
					targets_to_schedule[obs.name] = targetlist[priority == np.max(priority)][0]
					current_obs_tstart[i] = d
				else:
					# then other scopes with smaller FoV have to cover partial fields
					# might have to set aside second or third observations.  Kind of
					# need close to an integer ratio of 3 swope fields = 5 thacher fields
					# or something.  examining Tiling script will help.
					# build_tile_list() from Dave's tile.py script
					pass
					
		for obs in self.telescopes[np.argsort(fov)]:

			tel.write_schedule(obs,obs_date,targets_to_schedule[obs])


if __name__ == "__main__":

	usagestring = ""
	
	sched = Scheduler()

	parser = argparse.ArgumentParser(usage=usagestring, conflict_handler="resolve")
	parser.add_argument('-c','--configfile', default=None, type=str,
						help='configuration file')
	options, args = parser.parse_known_args()
		
	if options.configfile:
		config = configparser.ConfigParser()
		if not os.path.exists(options.configfile):
			raise RuntimeError('Configfile doesn\'t exist!')
		config.read(options.configfile)
	else:
		parser.print_help()
		raise RuntimeError('Configuration file must be specified at command line')

	parser = sched.add_options(usage=usagestring,config=config)
	args = parser.parse_args()

	sched.args = args
	sched.main()

