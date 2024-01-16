import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import gamma, riccati_jn, riccati_yn
from ChannelModelling.ChannelModel import calculate_kr_km_ka

print(calculate_kr_km_ka(2e6, 500e-9, 10e-6, complex(1.5, 0.0014)))

#print(miepython.mie(parameters["refraction"], 2*np.pi/parameters["wavelength"]*parameters["mean_radius"]))

"""
whp = []

whp_obj = cm.WHP([])
whp_obj.add_tuple(6, 2, 0.5)
whp_obj.add_tuple(6, 2, 0.5)
whp_bind = whp_obj.to_bind(0,1,3,5)
positions = [{"xt": 1, "yt": 0.5, "xr": 1, "yr": 1, "zi": 3, "zj": 1.8}]
dij = [np.sqrt((a["xt"]-a["xr"])**2+(a["yt"]-a["yr"])**2) for a in positions]
print(dij)
sh_param = cm.Shadowing_Parameters(positions, whp_bind)
sh_param_bind = sh_param.to_bind()
chb.calculate_expectancy(sh_param_bind)
"""