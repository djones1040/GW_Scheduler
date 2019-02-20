import configparser
import os

import Constants
from Observatory import Observatory
from Telescope import Swope, Nickel
from Utilities import *
from Target import TargetType, Target
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, ICRS, SkyOffsetFrame

from dateutil.parser import parse
import argparse
import astropy.coordinates as coord
from astropy import units as unit
from astropy import units as u

from datetime import datetime, timedelta
import time
from txtobj import txtobj
import numpy as np
from Tile import Tile
from shapely import geometry

class Scheduler:
	def __init__(self):
		self.telescopes = []

	def addwarning(self,warning):
		print(warning)
		self.warnings.append(warning)
		
	def add_options(self, parser=None, usage=None, config=None):
		if parser == None:
			parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

		parser.add_argument('-c','--configfile', default=None, type=str,
							help='configuration file')
		parser.add_argument("-f", "--file", default=config.get('main','file'),
							help="CSV file with targets to schedule.")
		parser.add_argument("-d", "--date", default=config.get('main','date'),
							help="YYYYMMDD formatted observation date.")
		parser.add_argument("--total_prob_goal", default=config.get('main','total_prob_goal'),type=float,
							help="total probability goal.  Helps to avoid low-prob targets with my simplified weighting scheme")
		parser.add_argument("--airmass_max", default=config.get('main','airmass_max'),type=float,
							help="maximum airmass to observe at")
		
		for telescope in config.sections():
			if telescope == 'main': continue
			self.telescopes += [telescope]
			
			parser.add_argument("--%s_slew_time_per_deg"%telescope,
								default=config.get(telescope,'slew_time_per_deg'),type=float)
			parser.add_argument("--%s_plate_scale"%telescope,
								default=config.get(telescope,'plate_scale'),type=float)
			parser.add_argument("--%s_detector_width_npix"%telescope,
								default=config.get(telescope,'detector_width_npix'),type=float)
			parser.add_argument("--%s_detector_height_npix"%telescope,
								default=config.get(telescope,'detector_height_npix'),type=float)
			parser.add_argument("--%s_lat"%telescope,
								default=config.get(telescope,'lat'),type=float)
			parser.add_argument("--%s_lon"%telescope,
								default=config.get(telescope,'lon'),type=float)
			parser.add_argument("--%s_elevation"%telescope,
								default=config.get(telescope,'elevation'),type=float)
			parser.add_argument("--%s_timezone"%telescope,
								default=config.get(telescope,'timezone'))
			parser.add_argument("--%s_gw_exptime"%telescope,
								default=config.get(telescope,'gw_exptime'),type=float)
			parser.add_argument("--%s_overhead_sec"%telescope,
								default=config.get(telescope,'overhead_sec'),type=float)
			parser.add_argument("--%s_obs_filter"%telescope,
								default=config.get(telescope,'obs_filter'))
			parser.add_argument("--%s_use"%telescope,default=config.get(telescope,'use'),
								help="schedule observations on this telescope?",type=bool)

		return parser

	def main(self):
	
		file_name = self.args.file
		obs_date = self.args.date

		observatories = {}
		targets_to_schedule = {}
		fov = []
		locations,aptelescopes = [],[]
		for t in self.telescopes:
			tel = Observatory(
				name=t,
				lon=self.args.__dict__['%s_lon'%t],
				lat=self.args.__dict__['%s_lat'%t],
				elevation=self.args.__dict__['%s_elevation'%t],
				horizon=-12,
				obs_date_str=obs_date,
				obs_filter=self.args.__dict__['%s_obs_filter'%t],
				timezone=self.args.__dict__['%s_timezone'%t],
				fov=self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]/60.,
				gw_exptime=self.args.__dict__['%s_gw_exptime'%t],
				slewtime_per_deg=self.args.__dict__['%s_slew_time_per_deg'%t]
			)
			tel.overhead = self.args.__dict__['%s_overhead_sec'%t]
			tel.detector_width = self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]
			tel.detector_height = self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]
			
			observatories[t] = tel
			targets_to_schedule[t] = []
			fov += [self.args.__dict__['%s_detector_width_npix'%t]*self.args.__dict__['%s_plate_scale'%t]/60.]

			location = EarthLocation.from_geodetic(tel.lon*u.deg, tel.lat*u.deg,
												   tel.elevation*u.m)
			locations += [location]
			aptelescopes += [Observer(location=location, timezone=tel.timezone)]

			
		target_data = txtobj(file_name,delimiter=',')

		names = target_data.__dict__['FieldID']
		ra = target_data.__dict__['R.A.[decimal]']
		dec = target_data.__dict__['Dec.[decimal]']
		if 'NetProb' in target_data.__dict__.keys():
			net_probs = target_data.__dict__['NetProb']
		else:
			net_probs = target_data.__dict__['2DProb']
		coords = SkyCoord(ra,dec,unit=(unit.hour, unit.deg))
		
		# sort by priorities
		iSort = np.argsort(net_probs)[::-1]
		names,ra,dec,net_probs = \
			names[iSort],ra[iSort],dec[iSort],net_probs[iSort]


		# datetime array
		def datetime_range(start, end, delta):
			current = start
			while current < end:
				yield current
				current += delta

		# TODO: make sure UTC offsets aren't so large that this screws up
		dts = [dt for dt in 
			   datetime_range(parse("%s 00:00" % obs_date),
							  parse("%s 00:00" % obs_date)+timedelta(hours=48), 
							  timedelta(seconds=5))]
		dts_airmass = [dt for dt in 
			datetime_range(parse("%s 00:00" % obs_date),
						   parse("%s 00:00" % obs_date)+timedelta(hours=48), 
						   timedelta(minutes=15))]

		
		probsum = 0
		for j in range(len(net_probs)):
			probsum += net_probs[j]
			if probsum > self.args.total_prob_goal: iProbMax = j

		targetdict = {}
		for t in self.telescopes:
			targetdict[t] = {}

		print('building target list')
		targets = np.array([])
		scheduled = np.zeros(len(names[0:iProbMax]))
		for n,c,p,l in zip(names[0:iProbMax],coords[0:iProbMax],net_probs[0:iProbMax],range(len(net_probs[0:iProbMax]))):

			target = Target(
				name=n, 
				coord=c, 
				net_prob=p,
				timezone=self.args.__dict__['%s_timezone'%t],
				observatory_lat=self.args.__dict__['%s_lat'%t],
				observatory_lon=self.args.__dict__['%s_lon'%t], 
				observatory_elev=self.args.__dict__['%s_elevation'%t], 
				sidereal_radian_array=observatories[t].sidereal_radian_array, 
				obs_date=observatories[t].obs_date,
				#obs_datetime=d,
				gw_exptime=self.args.__dict__['%s_gw_exptime'%t],
				obs_filter=self.args.__dict__['%s_obs_filter'%t],
				#tstart=d,
			)
			targets = np.append(targets,target)
			
			for t,l,apt in zip(self.telescopes,locations,aptelescopes):
				targetdict[t][n] = {}
				airmasslist = compute_airmass(
					c,apt,observatories[t].lon,observatories[t].lat,
					observatories[t].elevation,observatories[t].timezone,dts_airmass)
				targetdict[t][n]['airmass'] = list(np.concatenate([[am]*15*12 for am in airmasslist]))

				iRise = np.where((airmasslist > 0) & (airmasslist < 3))[0]
				if len(iRise): targetdict[t][n]['rise_time'] = dts_airmass[iRise[0]]
				else: targetdict[t][n]['rise_time'] = parse("2050/1/1")
				iSet = np.where(((airmasslist < 0) |
								 (airmasslist > 3)) &
								(np.array(dts_airmass) > targetdict[t][n]['rise_time']))[0]
				if len(iSet): targetdict[t][n]['set_time'] = dts[iSet[0]]
				else: targetdict[t][n]['set_time'] = parse("2050/1/1")

		print('finished building target list')
		
		telcoordlist = [None]*len(self.telescopes)
		current_obs_tstart = [None]*len(self.telescopes)
		current_obs_tend = [None]*len(self.telescopes)
		

		print('scheduling observations')
		for d,m in zip(dts,range(len(dts))):
			if np.sum(scheduled) == len(scheduled):
				continue
			for t,i in zip(
					np.array(self.telescopes)[np.argsort(fov)[::-1]],
					np.arange(len(self.telescopes))[np.argsort(fov)[::-1]]):

				# keep going if it's not night yet
				if observatories[t].utc_begin_night > d: continue
				if observatories[t].utc_end_night < d: continue

				if current_obs_tstart[i] and d < current_obs_tend[i]: #(d - current_obs_tstart[i]).seconds < observatories[t].gw_exptime:
					continue
				
				targetlist,coordlist,airmass_list,time_to_rise_list,time_to_set_list,slewtime_list,net_probs_list = \
					np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
				for target,n,c,p,l in zip(targets,names[0:iProbMax],coords[0:iProbMax],net_probs[0:iProbMax],range(len(net_probs[0:iProbMax]))):
				
					if scheduled[l]: continue
					airmass = targetdict[t][n]['airmass'][m]
					if airmass > 3 or airmass < 0: continue
					
					coordlist = np.append(coordlist,c)
					target.airmass = airmass
					net_probs_list = np.append(net_probs_list,p)
					
					if airmass < 3 and airmass >= 1:
						time_to_set = d - targetdict[t][n]['set_time']
						time_to_set_list = np.append(time_to_set_list,time_to_set.days*24*60*60 + time_to_set.seconds)
						if time_to_set.days > 1: target.time_to_set = 999
						else: target.time_to_set = time_to_set.days*24. + time_to_set.seconds/60./60.	
					else:
						time_to_set = 0
						target.time_to_set = 0
						time_to_set_list = np.append(time_to_set_list,time_to_set)
					
					if airmass > 3 or airmass < 1:
						time_to_rise = d - targetdict[t][n]['rise_time']
						time_to_rise_list = np.append(time_to_rise_list,time_to_rise.days*24*60*60 + time_to_rise.seconds)
						if time_to_rise.days > 1: target.time_to_rise = 999
						else: target.time_to_rise = time_to_rise.days*24. + time_to_rise.seconds/60./60.
					else:
						time_to_rise = 0
						target.time_to_rise = 0
						time_to_rise_list = np.append(time_to_rise_list,time_to_rise)

						
					
					if telcoordlist[i]:
						sep_deg = telcoordlist[i].separation(c).deg
						slewtime = observatories[t].slewtime_per_deg*sep_deg
					else:
						slewtime = 0

					targetlist = np.append(targetlist,target)
					airmass_list = np.append(airmass_list,airmass)
					slewtime_list = np.append(slewtime_list,slewtime)
					target.tstart = d + timedelta(seconds=slewtime)

					
				# compute modified priorities
				# if rising and AM > 2, lower priority (-1)
				# if setting w/i 30 min (AM > 2.5?), raise priority (+=3)
				# if setting w/i 60 min (AM > 2.5?), raise priority (+=2)
				# if slewtime < 15s, raise priority (+=1)
				# only look at tiles that comprise 50% of the total probability space
				# file should have a "done" flag, so you can get 50% of the total
				#	   probability w/o scheduling finished things.	Then you can
				#	   increase the threshold to 60%, 70%, etc as the night goes on

				if not len(targetlist):
					continue
				priority = np.ones(len(targetlist))
				priority[(airmass_list > 2) & (time_to_rise_list < 1)] -= 1
				priority[(time_to_set_list < 60)] += 2
				priority[(time_to_set_list < 30)] += 1
				priority[(slewtime_list < 15)] += 1
				priority += net_probs_list
				# then have to schedule things for x number of minutes
				
				if t == np.array(self.telescopes)[np.argsort(fov)[::-1]][0]:
					# TODO: FIX HACK!
					targetlist[priority == np.max(priority)][0].priority = np.max(priority)
					scheduled[names[0:iProbMax] == targetlist[priority == np.max(priority)][0].name] = 1
					targets_to_schedule[observatories[t].name] += [targetlist[priority == np.max(priority)][0]]
					current_obs_tstart[i] = d + timedelta(seconds=slewtime_list[priority == np.max(priority)][0]) + \
																   timedelta(seconds=observatories[t].overhead)
					current_obs_tend[i] = d + timedelta(seconds=slewtime_list[priority == np.max(priority)][0]) + \
																 timedelta(seconds=observatories[t].gw_exptime) + \
																 timedelta(seconds=observatories[t].overhead)

					telcoordlist[i] = coordlist[priority == np.max(priority)][0]
				else:
					# then other scopes with smaller FoV have to cover partial fields
					# might have to set aside second or third observations.	 Kind of
					# need close to an integer ratio of 3 swope fields = 5 thacher fields
					# or something.	 examining Tiling script will help.
					# build_tile_list() from Dave's tile.py script
					targetlist[priority == np.max(priority)][0].priority = np.max(priority)
					scheduled[names[0:iProbMax] == targetlist[priority == np.max(priority)][0].name] = 1

					telcoordlist[i] = coordlist[priority == np.max(priority)][0]

					targetcoordlist = build_tile_list_simple(
						telcoordlist[i],
						observatories[np.array(self.telescopes)[np.argsort(fov)[::-1]][0]].detector_width,
						observatories[np.array(self.telescopes)[np.argsort(fov)[::-1]][0]].detector_height,
						observatories[t].detector_width, observatories[t].detector_height)
					for targetcoord in targetcoordlist:
						target = Target(
							name=targetlist[priority == np.max(priority)][0].name,
							coord=targetcoord,
							net_prob=targetlist[priority == np.max(priority)][0].net_prob/len(targetcoordlist),
							timezone=self.args.__dict__['%s_timezone'%t],
							observatory_lat=self.args.__dict__['%s_lat'%t],
							observatory_lon=self.args.__dict__['%s_lon'%t], 
							observatory_elev=self.args.__dict__['%s_elevation'%t], 
							sidereal_radian_array=observatories[t].sidereal_radian_array, 
							obs_date=observatories[t].obs_date,
							gw_exptime=self.args.__dict__['%s_gw_exptime'%t],
							obs_filter=self.args.__dict__['%s_obs_filter'%t],
						)
						target.priority = np.max(priority)
						target.airmass = airmass
						target.time_to_rise = targetlist[priority == np.max(priority)][0].time_to_rise
						target.time_to_set = targetlist[priority == np.max(priority)][0].time_to_set
						
						targets_to_schedule[observatories[t].name] += [target]
					current_obs_tstart[i] = d + timedelta(seconds=slewtime_list[priority == np.max(priority)][0]) + \
																   timedelta(seconds=observatories[t].overhead)
					current_obs_tend[i] = d + timedelta(seconds=slewtime_list[priority == np.max(priority)][0]*(len(targetcoordlist)-1)) + \
																 timedelta(seconds=observatories[t].gw_exptime*len(targetcoordlist)) + \
																 timedelta(seconds=observatories[t].overhead*(len(targetcoordlist)-1))

					
		for t in np.array(self.telescopes)[np.argsort(fov)[::-1]]:
			write_schedule(targets_to_schedule[t],t)
		print('finished!')
			
def compute_airmass(coord, telescope, observatory_lon, observatory_lat,
					observatory_height, timezone,
					current_time):

		#target_coord = SkyCoord(coord.ra.deg,coord.dec.deg,unit=u.deg)
		airmass = telescope.altaz(current_time, coord).secz
		
		return airmass

			
def write_schedule(target_list,telescope_name):

	fout = open('%s_targets.txt'%telescope_name,'w')
	fout_sc = open('%s_targets.SkyCalc.txt'%telescope_name,'w')
	fout_ds9 = open('%s_targets.ds9.txt'%telescope_name,'w')
	print('# Object Name,Right Ascension,Declination,Filter,Exposure Time,Start Time,Net Prob,Priority,Airmass,HourstoRise,HourstoSet',file=fout)
	for target in target_list:
		print('%s,%s,%s,%s,%i,%s,%s,%s,%.2f,%.2f,%.2f'%(
			target.name,target.coord.ra.deg,target.coord.dec.deg,
			target.obs_filter,target.gw_exptime,target.tstart,
			target.net_prob,target.priority,target.airmass,
			target.time_to_rise,target.time_to_set),file=fout)
		print('%s %s %s 2000'%(
			target.name,GetSexigesimalString(target.coord.ra.deg,target.coord.dec.deg)[0],
			  GetSexigesimalString(target.coord.ra.deg,target.coord.dec.deg)[1]),file=fout_sc)
		print('%s %s %s 2000'%(
			target.name,target.coord.ra.deg,target.coord.dec.deg),file=fout_ds9)

	fout.close()
	fout_sc.close()
	fout_ds9.close()
	
def GetSexigesimalString(ra_decimal, dec_decimal):
    c = SkyCoord(ra_decimal,dec_decimal,unit=(u.deg, u.deg))
    ra = c.ra.hms
    dec = c.dec.dms

    ra_string = "%02d:%02d:%05.2f" % (ra[0],ra[1],ra[2])
    if dec[0] >= 0:
        dec_string = "+%02d:%02d:%05.2f" % (dec[0],np.abs(dec[1]),np.abs(dec[2]))
    else:
        dec_string = "%03d:%02d:%05.2f" % (dec[0],np.abs(dec[1]),np.abs(dec[2]))

    # Python has a -0.0 object. If the deg is this (because object lies < 60 min south), the string formatter will drop the negative sign
    if c.dec < 0.0 and dec[0] == 0.0:
        dec_string = "-00:%02d:%05.2f" % (np.abs(dec[1]),np.abs(dec[2]))
    return (ra_string, dec_string)

def build_tile_list_simple(cencoord, big_tile_width, big_tile_height,
						   tile_width, tile_height):

	dec_offset_half_big = coord.Angle(big_tile_height/2., unit=u.arcsec)
	ra_offset_half_big = coord.Angle(big_tile_width/2., unit=u.arcsec)

	mincoords = coord.SkyCoord(cencoord.ra.deg - ra_offset_half_big.degree,
							   cencoord.dec.deg - dec_offset_half_big.degree,
							   unit=(u.deg,u.deg))
	maxcoords = coord.SkyCoord(cencoord.ra.deg + ra_offset_half_big.degree,
							   cencoord.dec.deg + dec_offset_half_big.degree,
							   unit=(u.deg,u.deg))
		
	tilecoords = []
	
	n_tiles = (big_tile_width*big_tile_height)/(tile_width*tile_height)
	n_tiles_ra = big_tile_width/tile_width
	n_tiles_dec = big_tile_height/tile_height

	dec_offset_half = coord.Angle(tile_height/2., unit=u.arcsec)
	ra_offset_half = coord.Angle(tile_width/2., unit=u.arcsec)
	dec_offset = coord.Angle(tile_height, unit=u.arcsec)
	ra_offset = coord.Angle(tile_width, unit=u.arcsec)

	
	# TODO: rounding
	ra_tile = mincoords.ra.deg
	while ra_tile < maxcoords.ra.deg:
		dec_tile = mincoords.dec.deg
		while dec_tile < maxcoords.dec.deg:
			tilecoords += [coord.SkyCoord(ra_tile+ra_offset_half.degree,
										  dec_tile+dec_offset_half.degree,
										  unit=(u.deg,u.deg))]
			dec_tile = coord.SkyCoord(ra_tile,
									  dec_tile+dec_offset.degree,
									  unit=(u.deg,u.deg)).dec.deg

		ra_tile = coord.SkyCoord(ra_tile+ra_offset.degree,
								 dec_tile,
								 unit=(u.deg,u.deg)).ra.deg

	return tilecoords

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
