import pandas as pd
import numpy as np
import os, sys
from astropy import units as u
from astropy.coordinates import SkyCoord
import math

##################################### For the nearest velocity ########################################################

rcp_df = pd.read_csv('345.505+0.348rcp_imfit_pos.csv', sep='\t', usecols=[0,1,3,4,5], names=['vel', 'flux', 'Ra', 'Dec', 'offset'], header=None)
lcp_df = pd.read_csv('345.505+0.348lcp_imfit_pos.csv', sep='\t', usecols=[0,1,3,4,5], names=['vel', 'flux', 'Ra', 'Dec', 'offset'], header=None)

A_list_pairs = pd.read_csv('A_paired_sources.csv', usecols=[0,1,2,3], names=['Source Name', 'Trans', 'Vel_RCP', 'Vel_LCP'])


curr_path = os.getcwd().split('/')
tran_dir  = curr_path[-1][:4]
src_name  = curr_path[-1][5:]

src_name_list = A_list_pairs[(A_list_pairs['Source Name'] == src_name) & (A_list_pairs['Trans'] == tran_dir)]

rcp_Alist = src_name_list['Vel_RCP'].astype('float')
lcp_Alist = src_name_list['Vel_LCP'].astype('float')

rcp_imfit = rcp_df['vel'].astype('float')
lcp_imfit = lcp_df['vel'].astype('float')


def ang_sep(ra1, dec1, ra2, dec2):
	sep = np.round(math.acos(math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2)*math.cos(ra1-ra2))*3600, 1)

	return sep


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

final_rows = mul_frames.nsmallest(3, 'angular_offset')
print(final_rows)

final_rows.to_csv('zeeman_pairs.csv', encoding='utf-8')