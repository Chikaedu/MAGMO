import sys, os
import matplotlib.pyplot as plt
import matplotlib        as mpl
import matplotlib.gridspec as gridspec
import numpy as np
from numpy import loadtxt
from astropy.io import fits


data_P    = loadtxt('zeeman_positive_red.dat')
data_N    = loadtxt('zeeman_negative_blue.dat')
data_Null = loadtxt('non_split.dat')

Pos_x  = np.array(data_P[:, 0])
Pos_y  = np.array(data_P[:, 1])
Neg_x  = np.array(data_N[:, 0])
Neg_y  = np.array(data_N[:, 1])
Null_x = np.array(data_Null[:, 0])
Null_y = np.array(data_Null[:, 1])



###########################################################
###########################################################



data4 = loadtxt('spiralarm1.dat')
data5 = loadtxt('spiralarm2.dat')
data6 = loadtxt('spiralarm3.dat')
data7 = loadtxt('spiralarm4.dat')
data8 = loadtxt('spiralarm5.dat')


x4 = np.array(data4[:, 0])
y4 = np.array(data4[:, 1])

x5 = np.array(data5[:, 0])
y5 = np.array(data5[:, 1])

x6 = np.array(data6[:, 0])
y6 = np.array(data6[:, 1])

x7 = np.array(data7[:, 0])
y7 = np.array(data7[:, 1])

x8 = np.array(data8[:, 0])
y8 = np.array(data8[:, 1])






##########################################################
# end of data import
##########################################################


# PLOT

plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

mpl.rcParams['axes.linewidth'] = 2.


fig, ax = plt.subplots(2, sharex=True, figsize=(12,8))


# fig = plt.figure(figsize=(18,8))
# gs = gridspec.GridSpec(nrows=4, ncols=2)

# ax[0] = fig.add_subplot(gs[:2, 0])
# ax[1] = fig.add_subplot(gs[2:, 0])
# ax21 = fig.add_subplot(gs[:4, 1])
# fig.subplots_adjust(wspace=0.2, hspace=0.05)

# fts  = 22
# lbsz = 16
# lgds = 12


major_xticks = np.arange(-180., 60., 30.)
minor_xticks = np.arange(-180, 60., 15.)
# major_xticks = np.arange(60., -180., 30.)
# minor_xticks = np.arange(60, -180., 15.)

major_yticks = np.arange(-120, 100, 40.)
minor_yticks = np.arange(-120., 100., 20.)



## Axis 0
ax[0].plot(Null_x, Null_y, 'k*', markersize=4, markeredgewidth=1.5, alpha=0.4, fillstyle='none', label=r'$No$ $Splitting$')
ax[0].plot(Pos_x, Pos_y, 'ro', markersize=4, markeredgewidth=1.5, fillstyle='none', label=r'$Zeeman$ $positive$')
ax[0].plot(Neg_x, Neg_y, 'b+', markersize=4, markeredgewidth=1.5, fillstyle='none', label=r'$Zeeman$ $negative$')
ax[0].set_xlabel(r'\textbf{Longitude ($^{o}$)}', fontsize=16)
ax[0].set_ylabel(r'\textbf{Velocity (kms$^{-1}$)}', fontsize=16)
ax[0].set_xticks(major_xticks)
ax[0].set_xticks(minor_xticks, minor=True)
ax[0].set_yticks(major_yticks)
ax[0].set_yticks(minor_yticks, minor=True)
ax[0].tick_params(axis='x', labelsize=14, pad=8)
ax[0].tick_params(axis='y', labelsize=14)
ax[0].tick_params(which='both', width=1)
ax[0].tick_params(which='major', length=8, direction='in')
ax[0].tick_params(which='minor', length=4, direction='in')


ax[0].set_ylim(-150., 100.)
#ax[0].set_xlim(60, -180.)
ax[0].set_xlim(-180, 60.)
ax[0].legend(prop=dict(weight='bold', size=13.), loc='lower right')
ax[0].invert_xaxis()

# ax[0].legend(prop=dict(weight='bold', size=13), loc='lower right')



## Axis 1
ax[1].plot(Null_x, Null_y, 'k*', markersize=4, markeredgewidth=1.5, alpha=0.4, fillstyle='none', label=r'$No$ $Splitting$')
ax[1].plot(Pos_x, Pos_y,         'ro', markersize=4, markeredgewidth=1.5, fillstyle='none')
ax[1].plot(Neg_x, Neg_y,       'b+', markersize=4, markeredgewidth=1.5, fillstyle='none')
ax[1].plot(x4, y4,         'y-', lw=6, alpha=.45, label=r'$Norma$ $arm$')
ax[1].plot(x5, y5,         'b-', lw=6, alpha=.45, label=r'$Carina-Sagittarius$ $arm$')
ax[1].plot(x6, y6,         'r-', lw=6, alpha=.45, label=r'$Perseus$ $arm$')
ax[1].plot(x7, y7,         'g-', lw=6, alpha=.45, label=r'$Crux-Scutum$ $arm$')
ax[1].plot(x8, y8,         'c-', lw=6, alpha=.45, label=r'$Orion-Cygnus$')
ax[1].set_xlabel(r'\textbf{Longitude ($^{o}$)}', fontsize=16)
ax[1].set_ylabel(r'\textbf{Velocity (kms$^{-1}$)}', fontsize=16)
ax[1].set_xticks(major_xticks)
ax[1].set_xticks(minor_xticks, minor=True)
ax[1].set_yticks(major_yticks)
ax[1].set_yticks(minor_yticks, minor=True)
ax[1].tick_params(axis='x', labelsize=14, pad=8)
ax[1].tick_params(axis='y', labelsize=14)
ax[1].tick_params(which='both', width=1)
ax[1].tick_params(which='major', length=8, direction='in')
ax[1].tick_params(which='minor', length=4, direction='in')


ax[1].set_ylim(-150., 100.)
#ax[1].set_xlim(60, -180.)
ax[1].set_xlim(-180, 60.)
# ax[1].legend(prop=dict(weight='bold', size=13.), loc='lower right')
ax[1].invert_xaxis()




plt.savefig('Long_Vel.pdf', bbox_inches='tight', pad_inches=0.09, format='pdf')
plt.show()


















