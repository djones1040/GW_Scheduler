import Constants
from Observatory import Observatory
from Telescope import Swope, Nickel
from Utilities import *
from Target import TargetType, Target

from dateutil.parser import parse
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as unit

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file", help="CSV file with targets to schedule.")
	parser.add_argument("-d", "--date", help="YYYYMMDD formatted observation date.")
	parser.add_argument("-ot", "--obstele", help="Comma-delimited list of <Observatory>:<Telescope>, to schedule targets.")
	args = parser.parse_args()

	file_name = args.file
	obs_date = args.date
	observatory_telescopes = args.obstele.split(",")
	
	obs_keys = [o.split(":")[0] for o in observatory_telescopes]
	tele_keys = [t.split(":")[1] for t in observatory_telescopes]

	lco = Observatory(
		name="LCO",
		lon="-70.6915",
		lat="-29.0182",
		elevation=2402,
		horizon="-12",
		telescopes={"Swope":Swope()},
		obs_date_str=obs_date,
		utc_offset=lco_clt_utc_offset, # Chile observes Chile Standard Time (CLT) from 5/13/2017 - 8/12/2017 => UTC-4
		# utc_offset=lco_clst_utc_offset, # Chile observes Chile Summer Time (CLST) from 8/13/2017 - 12/31/2017 => UTC-3
		utc_offset_name="CLST"
	)

	lick = Observatory(
		name="Lick",
		lon="-121.6429",
		lat="37.3414",
		elevation=1283,
		horizon="-12",
		telescopes={"Nickel":Nickel()},
		obs_date_str=obs_date,
		utc_offset=lick_pdt_utc_offset, # California observes Pacific Daylight Time (PDT) from 3/12/2017 - 11/5/2017 => UTC-7
		# utc_offset=lick_pst_utc_offset, # California observes Pacific Standard Time (PST) from 1/1/2017 - 3/12/2017 => UTC-8
		utc_offset_name="PST"
	)

	observatories = {"LCO":lco, "Lick":lick}

	target_data = get_targets("%s" % file_name)
	names = [t[0] for t in target_data]
	ra = [t[1] for t in target_data]
	dec = [t[2] for t in target_data]
	priorities = [float(t[3]) for t in target_data]
	disc_dates = [t[4] for t in target_data]
	disc_mags = [float(t[5]) for t in target_data]
	types = [t[6] for t in target_data]
	coords = SkyCoord(ra,dec,unit=(unit.hour, unit.deg))

	for i in range(len(observatory_telescopes)):
		
		targets = []
		obs = observatories[obs_keys[i]]

		for j in range(len(names)):

			target_type = None
			disc_date = None

			if types[j] == "STD":
				target_type = TargetType.Standard
				disc_date = None
			elif types[j] == "TMP":
				target_type = TargetType.Template
				disc_date = parse(disc_dates[j])
			elif types[j] == "SN":
				target_type = TargetType.Supernova
				disc_date = parse(disc_dates[j])
			else:
				raise ValueError('Unrecognized target type!')

			targets.append(
				Target(
					name=names[j], 
					coord=coords[j], 
					priority=priorities[j], 
					target_type=target_type, 
					observatory_lat=obs.ephemeris.lat, 
					sidereal_radian_array=obs.sidereal_radian_array, 
					disc_date=disc_date, 
					apparent_mag=disc_mags[j], 
					obs_date=obs.obs_date
				)
			)

			obs.telescopes[tele_keys[i]].set_targets(targets)

		print("# of %s targets: %s" % (tele_keys[i], len(targets)))
		print("First %s target: %s" % (tele_keys[i], targets[0].name))
		print("Last %s target: %s" % (tele_keys[i], targets[-1].name))

		obs.schedule_targets(tele_keys[i])

	exit = input("\n\nENTER to exit")

if __name__ == "__main__": main()

		