import os, sys
import shutil
import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt
import matplotlib        as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import scipy as sp
from scipy import signal
import pickle
from uncertainties import ufloat
import re


parent_dir = os.path.basename(os.getcwd())
print(parent_dir)


src_dir = '1720'

output = open("Polarisation_percentage_list.dat", 'w')

#############################################################  Functions

def FindMax( listOfY , indexes, numOfMax):
    'This method finds the maximum values in the list of peaks'
    val = []
    listMax = []
    xList = []
    reconstructedList = []
    for c in range(0,numOfMax):

        if listOfY == val:
            listMax.append(0)
        else:
            listMax.append(max(listOfY))
            index = listOfY.index(max(listOfY))
            xList.append(indexes[index])
            listOfY.pop(index)
    return listMax, xList

def FindPeaks(listY):
    'This method finds the peaks from the list with the Y values'
    peaks = []
    indexes = []
    count = 0
    m2 = 0 #Old slope: starts with 0
    for value in range(1,len(listY)):
        m1 = listY[value] - listY[value-1] #New slope
        if( m2 > 0 and m1 < 0 ):
            peaks.append(listY[value-1])
            indexes.append( value-1 )
        m2 = m1 #Old slope is assigned
    return peaks, indexes


def calc_rms(amp): # returns rms
  x = len(amp)
  a = int(x / 10)
  rms_list = []
  for _set in range(9):
    rms = np.std(amp[(_set * a):(_set * a) + (2 * a)])
    rms_list.append(rms)
  median_rms = np.median(rms_list)
  return median_rms


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#############################################################


for file in sorted(os.listdir(src_dir)):
    src_path = os.path.join(src_dir, file)
    print(src_path)


    for maser in open('1665_median_velocity.txt').readlines():
        
        masers = maser.split()
        maser_name = str(masers[0])
        maser_vel = float(masers[3])

        if maser_name == file[:-10]:


            ##########################################################
            n = file[:-10]
            name = n.replace('.', '_')

            ############################################################


            data_q = loadtxt(src_path + '/' + 'Q.dat')
            data_u = loadtxt(src_path + '/' + 'U.dat')
            data_i = loadtxt(src_path + '/' + 'I.dat')
            data_v = loadtxt(src_path + '/' + 'V.dat')

            qy = np.array(data_q[:, 1])
            uy = np.array(data_u[:, 1])
            iy = np.array(data_i[:, 1])
            vy = np.array(data_v[:, 1])

            sig_i =  np.round(calc_rms(iy), 4)
            sig_v =  np.round(calc_rms(vy), 4)
            print('sigma I', sig_i)
            print('sigma V',sig_v)


            #############################################################
            #############################################################

            #    creating a fractional circular polarization file


            k = 0
            velocity = []
            stokes_i = []
            stokes_v = []
            stokes_q = []
            stokes_u = []

            output1 = open(src_path + '/' + "circ.dat", 'w')
            output2 = open(src_path + '/' + "err_circ.dat", 'w')
            output3 = open(src_path + '/' + "err_LIN.dat", 'w')


            for line in open(src_path + '/' + "I.dat").readlines():
                numbers = line.split()
                velocity.append(float(numbers[0]))
                stokes_i.append(float(numbers[1]))
                k = k+1

            for line in open(src_path + '/' + "V.dat").readlines():
                numbers = line.split()
                stokes_v.append(float(numbers[1]))

            for line in open(src_path + '/' + "Q.dat").readlines():
                numbers = line.split()
                stokes_q.append(float(numbers[1]))

            for line in open(src_path + '/' + "U.dat").readlines():
                numbers = line.split()
                stokes_u.append(float(numbers[1]))

            for j in range(0, k):
                circ = (stokes_v[j]/stokes_i[j])
                output1.write(str(velocity[j]) + '  ' + str(circ) + '\n')

                err_circ = np.round(np.abs(stokes_v[j])/np.abs(stokes_i[j]) *  np.sqrt( ((sig_v/np.abs(stokes_v[j]))**2)  + ((sig_i/np.abs(stokes_i[j]))**2)) * 100, 2)
                output2.write(str(velocity[j]) + '  ' + str(err_circ) + '\n')

                # err_circ = np.round( np.sqrt( ((sig_v/np.abs(stokes_i[j]))**2) ) * 100, 2)
                # output2.write(str(velocity[j]) + '  ' + str(err_circ) + '\n')

                err_linp = np.round((np.sqrt((np.abs(stokes_q[j])*np.abs(stokes_u[j]))/(np.abs(stokes_i[j]))**2) * 100), 2)
                output3.write(str(velocity[j]) + '  ' + str(err_linp) + '\n')


            output1.close()
            output2.close()
            output3.close()

            #############################################################
            #############################################################



            data = loadtxt(src_path + '/' + 'LIN.dat')
            data_c = loadtxt(src_path + '/' + 'circ.dat')

            linx = data[:, 0]
            liny = data[:, 1]
            xc = data_c[:, 0]
            yc = data_c[:, 1]



            sig_i =  np.round(calc_rms(iy), 4)
            sig_q =  np.round(calc_rms(qy), 4)
            sig_u =  np.round(calc_rms(uy), 4)
            sig_v =  np.round(calc_rms(vy), 4)
            sigma_qu = np.abs((sig_q + sig_u)/2)
            qu_err   = np.abs(sig_q*sig_u)

            sig_l = np.round(calc_rms(liny), 2)
            snr = sig_l*5
            list_less_snr = (liny < snr)
            liny[list_less_snr] = 0.


            ######################################################

            # making new files of circ and lin frac pol with signals greater than 5


            k = 0
            velocity = []
            stokes_i = []
            circ_pol = []
            stokes_lin = []


            output3 = open(src_path + '/' + "LIN_new.dat", 'w')

            for line in open(src_path + '/' + "I.dat").readlines():
                numbers = line.split()
                velocity.append(float(numbers[0]))
                stokes_i.append(float(numbers[1]))
                k = k+1

            for i in liny:
                stokes_lin.append(float(i))

            for j in range(0, k):

              output3.write(str(velocity[j]) + '  ' + str(stokes_lin[j]) + '\n')

            output3.close()


            ######################################################


            k = 0
            velocity = []
            stokes_i = []

            for line in open(src_path + '/' + "I.dat").readlines():
                numbers = line.split()
                velocity.append(float(numbers[0]))
                stokes_i.append(float(numbers[1]))
                k = k+1


            stokes_i = np.array(stokes_i)
            velocity = np.array(velocity)

            ul = find_nearest(velocity, maser_vel+10)
            ll = find_nearest(velocity, maser_vel-10)

            upper = int(np.where(velocity == ul)[0])
            lower = int(np.where(velocity == ll)[0])


            data_lin = loadtxt(src_path + '/' + 'LIN_new.dat')
            xl = data[:, 0]
            yl = data[:, 1]

            ################################################################ using stokes I to mask the circ and lin frac pols
            stokes_i = stokes_i[lower:upper] 
            indx    = np.where( stokes_i < 5. * sig_i )[0]
            linp           = yl[lower:upper]

            if max(linp) == 0.:

                linp_nan           = linp + 2.*sigma_qu
                LIN_l_4        = (linp_nan / stokes_i)*100.
                LIN_l_4[indx]  = 0.

            elif max(linp) > 0.:

                linp           = ((linp**2 - sigma_qu**2))**0.5
                linp_nan       = np.nan_to_num(linp)
                LIN_l_4        = (linp_nan / stokes_i)*100.
                LIN_l_4[indx]  = 0.
            else:

                continue

            # print('final: ', LIN_l_4)

            stokes_i[indx] = 0.
            # ################################################################
            circ          = np.abs(yc)
            circp         = circ[lower:upper]
            circ_4        = (circp *100.)
            circ_4[indx]  = 0.
            ################################################################

            velocity = velocity[lower:upper]


            # make an example object to pickle
            xobj = {'x':velocity,
                    'yl':LIN_l_4,
                    'yc':circ_4,
                    'foo':True,
                    'spam':False}
            pickle.dump(xobj, open(src_path + '/' + 'polar_plot.pickle', 'wb'))
            # Load data
            dat = pickle.load(open(src_path + '/' + 'polar_plot.pickle', "rb"))
            x = dat['x']
            yl = dat['yl']
            yc = dat['yc']


            # Plot 


            plt.rc('font', weight='bold')
            plt.rc('text', usetex=True)
            plt.rc('xtick', labelsize=12)
            plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
            mpl.rcParams['axes.linewidth'] = 1.5
            legend_properties = {'weight':'bold'}
            mpl.rcParams["legend.borderaxespad"] = 0.8



            fig, ax1 = plt.subplots()

            lns1 = ax1.plot(x, yc, 'r-', label=r'\textbf{Circular}', lw=1.5)
            lns2 = ax1.plot(x, yl, 'b--', label=r'\textbf{Linear}', lw=1.5)
            plt.tight_layout()
            ax1.set_xlabel(r'\textbf{LSR Velocity (km $s^{-1}$)}', fontsize=12)
            ax1.set_ylabel(r'\textbf{Fractional polarisation (\%)}', fontsize=14)
            ax1.tick_params(axis='x', labelsize=14)
            ax1.tick_params(axis='y', labelsize=14)
            ax1.tick_params(which='both', width=2)
            ax1.tick_params(which='major', length=4, direction='in')
            ax1.set_xlim(maser_vel-10, maser_vel+10)
            # ax1.set_ylim(ymin=-10.)
            # ax1.yaxis.set_major_locator(MultipleLocator(20))
            plt.title(r'\textbf{' + n + ' ' + '(' + src_dir + '-MHz' + ')' +'}')


            ax2      = ax1.twinx()

            lns3 = ax2.plot(x, stokes_i, 'k-', label=r'\textbf{Stokes-I}', lw=4.0, alpha=0.5)
            ax2.set_ylabel(r'\textbf{Total flux density (Jy)}', fontsize=14)
            align_yaxis(ax1, 0, ax2, 0)
            # ax2.set_ylim(ymin=0.)
            ax2.tick_params(axis='y', labelsize=14)
            ax2.tick_params(which='both', width=2)
            ax2.tick_params(which='major', length=4, direction='in')


            lns = lns1+lns2+lns3
            labs = [l.get_label() for l in lns]
            ax2.legend(lns, labs, loc=0)


            plt.savefig(src_path + '/' + n + 'percent_frac_pols.pdf', format = 'pdf', dpi=plt.gcf().dpi, bbox_inches = 'tight')
            plt.savefig('polarisation_' + src_dir + '/' + name + '_' + src_dir + '.pdf', format = 'pdf', dpi=plt.gcf().dpi, bbox_inches = 'tight')
            # plt.show()
            plt.close()




            ########################################################  maximum of circ and lin pol

            peaksList_circ = FindPeaks(yc)
            peaksList_lin  = FindPeaks(yl)


            # print(peaksList_circ)
            # print(peaksList_lin)


            data_err_circ = loadtxt(src_path + '/' + 'err_circ.dat')
            data_err_lin  = loadtxt(src_path + '/' + 'err_LIN.dat')
            err_c         = data_err_circ[:, 1]
            err_l         = data_err_lin[:, 1]



            iy    = iy[lower:upper]
            vy    = vy[lower:upper]
            qy    = qy[lower:upper]
            uy    = uy[lower:upper]
            liny  = liny[lower:upper]
            err_c = err_c[lower:upper]
            err_l = err_l[lower:upper]
            # print('heres length: ', len(yc))
            # print('heres max: ',max(yc))



            if peaksList_circ and peaksList_lin[0] == []:

                print('only lin empty')

                maxList_circ = FindMax(peaksList_circ[0],peaksList_circ[1],1)
                circular = float(str(maxList_circ[0])[1:-1])
                max_circ = np.round(circular)
                vel_index_circ =  int(str(maxList_circ[1])[1:-1])

                err_circ = round(err_c[vel_index_circ])

                circular   = ufloat(max_circ, err_circ)
                output.write(n + '   ' + str("{:.1f}".format(circular)) + '   ' + str('<5') + "\n")



            elif peaksList_circ and peaksList_lin:

                print('both not empty')


                maxList_circ = FindMax(peaksList_circ[0],peaksList_circ[1],1)
                maxList_lin = FindMax(peaksList_lin[0],peaksList_lin[1],1)
                circular = float(str(maxList_circ[0])[1:-1])
                linear   = float(str(maxList_lin[0])[1:-1])
                max_circ = np.round(circular)
                max_lin  = np.round(linear)

                vel_index_circ =  int(str(maxList_circ[1])[1:-1])
                vel_index_lin  =  int(str(maxList_lin[1])[1:-1])

                err_circ     = round(err_c[vel_index_circ])
                geom_err_lin = round(err_l[vel_index_circ])

                circular   = ufloat(max_circ, err_circ)
                linear     = ufloat(max_lin, geom_err_lin)


                output.write(n + '   ' + str("{:.1f}".format(circular)) + '   ' + str("{:.1f}".format(linear)) + "\n")


            elif peaksList_circ[0] == [] and peaksList_lin[0] == []:

                print('both are empty')

                output.write(n + '   ' + str('<5') + '   ' + str('<5') + "\n")

            else:

                pass

output.close()
sys.exit()

