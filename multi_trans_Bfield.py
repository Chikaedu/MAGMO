import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy as np
from numpy import loadtxt
import scipy.stats as stats
import pandas as pd


data = loadtxt('Desktop/lv_data/multi_trans.dat')

tot_pair = np.array(data[:, 0])
eff_dir = np.array(data[:, 1])
normed = np.array(data[:, 2])

# print(len(normed))

# pos_one = (np.abs(normed) < 1.0).all() and (np.abs(normed) > 0.0).all()
# print(len(normed[pos_one]))


bins = 7

# plt.hist(normed, bins)
# # plt.hist(np.abs(normed), bins)
# # plt.xlim(-1, 1)
# plt.show()


lower, upper = -1.0, 1.0
sample_size = 78
mean = 0.0
std_dev = 0.5

# random_sample = np.random.normal(loc=mean, scale=std_dev, size=sample_size)


boundary = stats.truncnorm(
    (lower - mean) / std_dev, (upper - mean) / std_dev, loc=mean, scale=std_dev)

random_sample = stats.norm(loc=mean, scale=std_dev)


# fig, ax = plt.subplots(2, sharex=False)
# ax[0].hist(boundary.rvs(78), bins, histtype='step', color='k', label='fake', density=False)
# ax[0].legend()
# ax[1].hist(normed, bins, histtype='step', color='r', label='real', density=False)
# ax[1].legend()
# plt.savefig('probability_distribution_2', format='pdf')
# plt.show()


fig, ax = plt.subplots(1, sharex=False)
ax.hist(np.abs(boundary.rvs(78)), bins, histtype='step', color='k', label='fake', density=False)
ax.hist(np.abs(normed), bins, histtype='step', color='r', label='real', density=False)
ax.legend(loc='upper left')
# ax[1].hist(normed, bins, histtype='step', color='r', label='real', density=False)
# ax[1].legend()
plt.savefig('probability_distribution_b', format='pdf')
plt.show()
