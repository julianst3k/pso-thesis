import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import logging
from numpy import genfromtxt
from ChannelSimulation.plots.plot_utils import open_file, do_exponential_form, H_wrapper
x_max = 5
x_min = 0
y_max = 3
y_min = 1
height = 1.8
time = 1
intensity = 5
arr = open_file('../at_{}_user_{}_normd.csv'.format(1, height), (100,100))
arr_exp = do_exponential_form(arr, intensity, time)
arr2 = open_file('../at_{}_user_{}_normd.csv'.format(2, height), (100,100))
np_x = np.linspace(x_min, x_max, 100)
np_y = np.linspace(y_min, y_max, 100)
zr = 1.8
zt = 3
xr = 1.5
yr = 2
yt = 0.5
xt = 1


X, Y = np.meshgrid(np_x, np_y)
pab, pap, me = H_wrapper(X, Y, xt, yt, xr, yr, zr, zt)

fig, ax = plt.subplots(1,3,subplot_kw={"projection": "3d"})
surf = ax[0].plot_surface(X, Y, pab, cmap = cm.coolwarm)
fig.colorbar(surf, shrink = 0.25, aspect = 10, location= "left")

surf = ax[1].plot_surface(X, Y, pap, cmap = cm.coolwarm)
fig.colorbar(surf, shrink = 0.25, aspect = 10, location= "left")

surf = ax[2].plot_surface(X, Y, me, cmap = cm.coolwarm)
fig.colorbar(surf, shrink = 0.25, aspect = 10, location= "left")
for i in range(3):
    ax[i].set_xlabel("X [m]")
    ax[i].set_ylabel("Y [m]")
    ax[i].set_title("Minimum Height for Obstruction.")
ax[0].set_zlabel("Minimum Height [m]")
plt.show()