import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import logging
from numpy import genfromtxt
from ChannelSimulation.plots.plot_utils import open_file, do_exponential_form, H_wrapper

arr = open_file('../response_file.csv', (8,4,4,180))

fig, ax = plt.subplots(1,4)
for i in range(4):
    ax[i].plot(arr[1,i,0,:],)
    ax[i].plot(arr[1,i,1,:])
    ax[i].plot(arr[1,i,2,:])
    ax[i].plot(arr[1,i,3,:])

plt.show()



