import os, sys
import shutil
from mirpy import miriad
import numpy as np
from numpy import loadtxt
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import math


def uvlist_filt(output):
    return output.split('\n')

miriad.set_filter('imfit',uvlist_filt)
miriad.set_filter('maxfit',uvlist_filt)


curr_path = os.getcwd().split('/')
tran_dir  = curr_path[-1][:4]


src_dir = "source1/"

for file in os.listdir(src_dir):


	A_list_pairs = pd.read_csv('~/Desktop/lv_data/A_paired_sources.csv', usecols=[0,1,2,3], names=['Source Name', 'Trans', 'Vel_RCP', 'Vel_LCP'])
	src_name_list = A_list_pairs[(A_list_pairs['Source Name'] == file) & (A_list_pairs['Trans'] == tran_dir)]

	rcp_Alist = src_name_list['Vel_RCP'].astype('float')
	lcp_Alist = src_name_list['Vel_LCP'].astype('float')	
	vel_range_list = rcp_Alist.tolist() + lcp_Alist.tolist()

	OH_min_vel, OH_max_vel = min(vel_range_list)-2, min(vel_range_list)+2            #  velocity range of maser source. This is obtained from Table 1 of MAGMO paper 2.
	print(OH_min_vel, OH_max_vel)
	#							 													 #	WIll be used to obtain positions at certain velocity intervals of the detection. 
	chan_num = 1																	 #  Reference or pivot channel number of the image cube, in this case I use the first channel.


	src_path = os.path.join(src_dir, file)

	vel_chan_1  = round(float(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor', region = 'percentage(1,1)('+str(chan_num)+')')[15][28:41]), 2)
	Ref_Ra   	= 		str(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor', region = 'box(256,256,256,256)(1)')[13][28:])
	Ref_Dec  	= 		str(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor', region = 'box(256,256,256,256)(1)')[14][28:])
	Ref_coord   =       SkyCoord(Ref_Ra, Ref_Dec, unit=(u.hourangle, u.deg), frame='icrs')


	print(vel_chan_1)

	if tran_dir == ('1665') or tran_dir == ('1665'):
		min_chan  = int((np.abs(vel_chan_1 - OH_min_vel)/0.088))
		max_chan  = int((np.abs(vel_chan_1 - OH_max_vel)/0.088))
	elif tran_dir == '1612':
		min_chan  = int((np.abs(vel_chan_1 - OH_min_vel)/0.091))
		max_chan  = int((np.abs(vel_chan_1 - OH_max_vel)/0.091))
	elif tran_dir == '1720':
		min_chan  = int((np.abs(vel_chan_1 - OH_min_vel)/0.085))
		max_chan  = int((np.abs(vel_chan_1 - OH_max_vel)/0.085))
	else:
		pass

	print(min_chan, max_chan)

	output_rcp = open(src_path + '/' + file + 'rcp_imfit_pos.csv', 'w')

	for chan in range(min_chan, max_chan):
			Src_Ra   	 = 		str(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor', region = 'percentage(2,2)(' + str(chan) + ')')[13][28:])
			Src_Dec  	 = 		str(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor', region = 'percentage(2,2)(' + str(chan) + ')')[14][28:])
			max_vel      = round(float(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor',  region = 'percentage(2,2)(' + str(chan) + ')')[15][28:41]), 2)
			max_flux     = round(float(miriad.maxfit( _in = src_path + '/' + file + 'rcp_irestor',  region = 'percentage(2,2)(' + str(chan) + ')')[3][31:]), 2)
			OH_coord     = SkyCoord(Src_Ra, Src_Dec, unit=(u.hourangle, u.deg), frame='icrs')
			OH_name      = OH_coord.galactic
			ang_sep      = round(OH_coord.separation(Ref_coord).arcsecond, 1)
			lon          = round(OH_name.l.degree, 6)
			lat          = round(OH_name.b.degree, 6)

			if lat < 0:
				output_rcp.write(str(max_vel) + '	' + str(max_flux) + '	' + file + '	' + str(lon) + '	' + str(lat) + '	' + str(ang_sep) + '   ' + Src_Ra + '	' + Src_Dec + "\n")
			else:
				output_rcp.write(str(max_vel) + '	' + str(max_flux) + '	' + file + '	' + str(lon) + '	' + str(lat) + '	' + str(ang_sep) + '	' + Src_Ra + '	' + Src_Dec + "\n")
	output_rcp.close()


	output_lcp = open(src_path + '/' + file + 'lcp_imfit_pos.csv', 'w')

	for chan in range(min_chan, max_chan):

			Src_Ra   	 = 		str(miriad.maxfit( _in = src_path + '/' + file + 'lcp_irestor', region = 'percentage(2,2)(' + str(chan) + ')')[13][28:])
			Src_Dec  	 = 		str(miriad.maxfit( _in = src_path + '/' + file + 'lcp_irestor', region = 'percentage(2,2)(' + str(chan) + ')')[14][28:])
			max_vel      = round(float(miriad.maxfit( _in = src_path + '/' + file + 'lcp_irestor',  region = 'percentage(2,2)(' + str(chan) + ')')[15][28:41]), 2)
			max_flux     = round(float(miriad.maxfit( _in = src_path + '/' + file + 'lcp_irestor',  region = 'percentage(2,2)(' + str(chan) + ')')[3][31:]), 2)
			OH_coord     = SkyCoord(Src_Ra, Src_Dec, unit=(u.hourangle, u.deg), frame='icrs')
			OH_name      = OH_coord.galactic
			ang_sep      = round(OH_coord.separation(Ref_coord).arcsecond, 1)
			lon          = round(OH_name.l.degree, 6)
			lat          = round(OH_name.b.degree, 6)

			if lat < 0:
				output_lcp.write(str(max_vel) + '	' + str(max_flux) + '	' + file + '	' + str(lon) + '	' + str(lat) + '	' + str(ang_sep) + '   ' + Src_Ra + '	' + Src_Dec + "\n")
			else:
				output_lcp.write(str(max_vel) + '	' + str(max_flux) + '	' + file + '	' + str(lon) + '	' + str(lat) + '	' + str(ang_sep) + '	' + Src_Ra + '	' + Src_Dec + "\n")
	output_lcp.close()


	##################################### For the nearest velocity ########################################################

	rcp_df = pd.read_csv(src_path + '/' + file + 'rcp_imfit_pos.csv', sep='\t', usecols=[0,1,3,4,5], names=['vel', 'flux', 'Ra', 'Dec', 'offset'], header=None)
	lcp_df = pd.read_csv(src_path + '/' + file + 'lcp_imfit_pos.csv', sep='\t', usecols=[0,1,3,4,5], names=['vel', 'flux', 'Ra', 'Dec', 'offset'], header=None)
	A_list_pairs = pd.read_csv('~/Desktop/lv_data/A_paired_sources.csv', usecols=[0,1,2,3], names=['Source Name', 'Trans', 'Vel_RCP', 'Vel_LCP'])

	src_name_list = A_list_pairs[(A_list_pairs['Source Name'] == file) & (A_list_pairs['Trans'] == tran_dir)]
	rcp_Alist = src_name_list['Vel_RCP'].astype('float')
	lcp_Alist = src_name_list['Vel_LCP'].astype('float')
	rcp_imfit = rcp_df['vel'].astype('float')
	lcp_imfit = lcp_df['vel'].astype('float')


	mid_ra_rcp  = []
	mid_dec_rcp = []
	mid_ra_lcp  = []
	mid_dec_lcp = []
	up_ra_rc    = []
	up_dec_rc   = []
	up_ra_lc    = []
	up_dec_lc   = []
	low_ra_rc   = []
	low_dec_rc  = []
	low_ra_lc   = []
	low_dec_lc  = []


	################### Middle row ######################

	list_middle_rc = list()
	for i in rcp_Alist:
		ind_vel = min(rcp_imfit, key=lambda x:abs(x-(i)))
		mid_ra_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][3])
		list_middle_rc.append(rcp_df[rcp_df['vel'] == ind_vel])
	middle_df_rc = pd.concat(list_middle_rc)
	middle_df_rc.columns = ['vel_rc', 'flux_rc', 'Ra_rc', 'Dec_rc', 'offset_rc']

	list_middle_lc = list()
	for i in lcp_Alist:
		ind_vel = min(lcp_imfit, key=lambda x:abs(x-(i)))
		mid_ra_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][3])
		list_middle_lc.append(lcp_df[lcp_df['vel'] == ind_vel])
	middle_df_lc = pd.concat(list_middle_lc)
	middle_df_lc.columns = ['vel_lc', 'flux_lc', 'Ra_lc', 'Dec_lc', 'offset_lc']

	middle_df_rc.reset_index(drop=True, inplace=True)
	middle_df_lc.reset_index(drop=True, inplace=True)

	middle_df = pd.concat([middle_df_rc,middle_df_lc], axis=1)
	pd.set_option("display.max_rows", None, "display.max_columns", None)


	mid_ang_off = np.sin(middle_df['Dec_rc'])*np.sin(middle_df['Dec_lc']) + np.cos(middle_df['Dec_rc'])*np.cos(middle_df['Dec_lc'])*np.cos(middle_df['Ra_rc']-middle_df['Ra_lc'])
	middle_df['angular_offset'] = np.round(np.arccos(mid_ang_off)*3600, 1)


	################### Upper row ######################
	list_upper_rc = list()
	for i in rcp_Alist:
		ind_vel = min(rcp_imfit, key=lambda x:abs(x-(i-0.09)))
		mid_ra_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][3])
		list_upper_rc.append(rcp_df[rcp_df['vel'] == ind_vel])
	upper_df_rc = pd.concat(list_upper_rc)
	upper_df_rc.columns = ['vel_rc', 'flux_rc', 'Ra_rc', 'Dec_rc', 'offset_rc']

	list_upper_lc = list()
	for i in lcp_Alist:
		ind_vel = min(lcp_imfit, key=lambda x:abs(x-(i-0.09)))
		mid_ra_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][3])
		list_upper_lc.append(lcp_df[lcp_df['vel'] == ind_vel])
	upper_df_lc = pd.concat(list_upper_lc)
	upper_df_lc.columns = ['vel_lc', 'flux_lc', 'Ra_lc', 'Dec_lc', 'offset_lc']

	upper_df_rc.reset_index(drop=True, inplace=True)
	upper_df_lc.reset_index(drop=True, inplace=True)

	upper_df = pd.concat([upper_df_rc,upper_df_lc], axis=1)
	pd.set_option("display.max_rows", None, "display.max_columns", None)


	upper_ang_off = np.sin(upper_df['Dec_rc'])*np.sin(upper_df['Dec_lc']) + np.cos(upper_df['Dec_rc'])*np.cos(upper_df['Dec_lc'])*np.cos(upper_df['Ra_rc']-upper_df['Ra_lc'])
	upper_df['angular_offset'] = np.round(np.arccos(upper_ang_off)*3600, 1)


	################### Lower row ######################

	list_lower_rc = list()
	for i in rcp_Alist:
		ind_vel = min(rcp_imfit, key=lambda x:abs(x-(i+0.09)))
		mid_ra_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_rcp.append(rcp_df.loc[rcp_df['vel'] == ind_vel].iloc[0][3])
		list_lower_rc.append(rcp_df[rcp_df['vel'] == ind_vel])
	lower_df_rc = pd.concat(list_lower_rc)
	lower_df_rc.columns = ['vel_rc', 'flux_rc', 'Ra_rc', 'Dec_rc', 'offset_rc']

	list_lower_lc = list()
	for i in lcp_Alist:
		ind_vel = min(lcp_imfit, key=lambda x:abs(x-(i+0.09)))
		mid_ra_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][2])
		mid_dec_lcp.append(lcp_df.loc[lcp_df['vel'] == ind_vel].iloc[0][3])
		list_lower_lc.append(lcp_df[lcp_df['vel'] == ind_vel])
	lower_df_lc = pd.concat(list_lower_lc)
	lower_df_lc.columns = ['vel_lc', 'flux_lc', 'Ra_lc', 'Dec_lc', 'offset_lc']

	lower_df_rc.reset_index(drop=True, inplace=True)
	lower_df_lc.reset_index(drop=True, inplace=True)

	lower_df = pd.concat([lower_df_rc,lower_df_lc], axis=1)
	pd.set_option("display.max_rows", None, "display.max_columns", None)


	lower_ang_off = np.sin(lower_df['Dec_rc'])*np.sin(lower_df['Dec_lc']) + np.cos(lower_df['Dec_rc'])*np.cos(lower_df['Dec_lc'])*np.cos(lower_df['Ra_rc']-lower_df['Ra_lc'])
	lower_df['angular_offset'] = np.round(np.arccos(lower_ang_off)*3600, 1)


	########################################################################

	frames = [middle_df, upper_df, lower_df]
	mul_frames = pd.concat(frames, keys=['mid','upp','low'])
	print(mul_frames)

	final_rows = mul_frames.nsmallest(len(rcp_Alist), 'angular_offset')
	print(final_rows)

	final_rows.to_csv(src_path + '/' + file + 'zeeman_pairs.csv', encoding='utf-8')


sys.exit()
