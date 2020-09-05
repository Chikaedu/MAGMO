import os, sys
from mirpy import miriad
import numpy as np
from numpy import loadtxt
import shutil
from astropy import units as u
from astropy.coordinates import SkyCoord



def uvlist_filt(output):
    return output.split('\n')

miriad.set_filter('imfit',uvlist_filt)
miriad.set_filter('maxfit',uvlist_filt)


pivot_Src_Ra, pivot_Src_Dec  = '16h09m52.37s', '-51d54m57.60s'
OH_min_vel, OH_max_vel       = -100, -80


src_dir = "1667/"

for file in os.listdir(src_dir):

	src_path = os.path.join(src_dir, file)

	for meth in open('MMB_list_com.txt').readlines():

		meths     = meth.split()
		meth_name = str(meths[0])
		meth_Ra   = str(meths[2])
		meth_Dec  = str(meths[3])

		if meth_name == file:

			output_properties = open(src_path + '/' + file + 'mutiple_positions', 'w')

			mid_chan_num = 852
			max_vel  = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[15][28:41]), 2)

			print(max_vel)

			min_chan_off, max_chan_off  = int(round((max_vel-OH_min_vel)/0.088)), int(round((OH_max_vel-max_vel)/0.088))
			min_chan_num, max_chan_num  = mid_chan_num-min_chan_off, mid_chan_num+max_chan_off
			tot_chan_num = min_chan_off+max_chan_off
			print(min_chan_num, max_chan_num)
			print(tot_chan_num)


			for chan in range(min_chan_num, max_chan_num):


				Src_Ra   	= 		str(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = 'percentage(5,5)(' + str(chan) + ')')[13][28:])
				Src_Dec  	= 		str(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = 'percentage(5,5)(' + str(chan) + ')')[14][28:])
				max_vel  = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr',  region = 'percentage(5,5)(' + str(chan) + ')')[15][28:41]), 2)
				max_flux = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr',  region = 'percentage(5,5)(' + str(chan) + ')')[3][31:]), 2)

				OH_coord        = SkyCoord(Src_Ra, Src_Dec, unit=(u.hourangle, u.deg), frame='icrs')
				pivot_coord     = SkyCoord(pivot_Src_Ra, pivot_Src_Dec, unit=(u.hourangle, u.deg), frame='icrs')
				OH_name         = OH_coord.galactic
				meth_coord      = SkyCoord(meth_Ra, meth_Dec, unit=(u.hourangle, u.deg), frame='icrs')
				meth_ang_sep    = round(OH_coord.separation(meth_coord).arcsecond, 1)
				OH_ang_sep      = round(OH_coord.separation(pivot_coord).arcsecond, 1)
				lon             = round(OH_name.l.degree, 3)
				lat             = round(OH_name.b.degree, 3)


				if lat < 0:

					output_properties.write(file + '	' + str(lon) +  str(lat) + '	' + str(meth_ang_sep) + '	' + str(OH_ang_sep) + '	' + Src_Ra + '	' + Src_Dec + '	' + str(max_vel) + '	' + str(max_flux) + "\n")

				else:
					output_properties.write(file + '	' + str(lon) +  str(lat) + '	' + str(meth_ang_sep) + '	' + str(OH_ang_sep) + '	' + Src_Ra + '	' + Src_Dec + '	' + str(max_vel) + '	' + str(max_flux) + "\n")

			output_properties.close()
sys.exit()














































# src_dir = "source1/"

# for file in os.listdir(src_dir):

# 	src_path = os.path.join(src_dir, file)

# 	for meth in open('MMB_list_com.txt').readlines():

# 		meths     = meth.split()
# 		meth_name = str(meths[0])
# 		meth_vel  = float(meths[1])
# 		meth_Ra   = str(meths[2])
# 		meth_Dec  = str(meths[3])

# 		if meth_name == file:


# 									## maxfit

# 			Src_Ra   	= 		str(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[13][28:])
# 			Src_Dec  	= 		str(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[14][28:])
# 			Src_Amp  = round(float(miriad.imfit(  _in = src_path + '/' + file + '.irestor', object = 'point', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[15][30:41]), 2)
# 			Amp_err  = round(float(miriad.imfit(  _in = src_path + '/' + file + '.irestor', object = 'point', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[15][46:]), 2)
# 			Beam_Pa  =       float(miriad.imfit(  _in = src_path + '/' + file + '.irestor', object = 'point', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[9][38:44])
# 			Bm_Maj   =             miriad.imfit(  _in = src_path + '/' + file + '.irestor', object = 'point', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[8][37:45]
# 			Bm_Min   =             miriad.imfit(  _in = src_path + '/' + file + '.irestor', object = 'point', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[8][48:54]
# 			RA_offset   = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[8][31:47]), 2)
# 			Dec_offset  = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[9][31:47]), 2)
# 			max_vel  = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[15][28:41]), 2)
# 			max_flux = round(float(miriad.maxfit( _in = src_path + '/' + file + '.ipbcorr', region = region = 'percentage(5,5)(' + str(mid_chan_num) + ')')[3][30:]), 2)



# 			OH_coord   = SkyCoord(Src_Ra, Src_Dec, unit=(u.hourangle, u.deg), frame='icrs')
# 			OH_name    = OH_coord.galactic
# 			meth_coord = SkyCoord(meth_Ra, meth_Dec, unit=(u.hourangle, u.deg), frame='icrs')
# 			ang_sep    = round(OH_coord.separation(meth_coord).arcsecond, 1)
# 			lon        = round(OH_name.l.degree, 3)
# 			lat        = round(OH_name.b.degree, 3)



# 			if lat < 0:

# 				output_properties.write(file + '	' + str(lon) +  str(lat) + '	' + str(ang_sep) + '	' + Src_Ra + '	' + Src_Dec + '	' + str(max_vel) + '	' + str(max_flux))

# 			else:
# 				output_properties.write(file + '	' + str(lon) +  str(lat) + '	' + str(ang_sep) + '	' + Src_Ra + '	' + Src_Dec + '	' + str(max_vel) + '	' + str(max_flux))

# 				output_properties.close()

# sys.exit()
