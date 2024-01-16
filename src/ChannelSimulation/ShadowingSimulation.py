import ChannelModelling.ChannelModel as cm
import model as chb
import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import logging

FORMAT = '%(asctime)s %(name)-12s %(message)s'
logging.basicConfig(format=FORMAT,level=logging.INFO,force=True)
logger = logging.getLogger('Shadowing')

class WHP_Wrapper:
    def __init__(self):
        self.whp_obj = cm.WHP([])

    def add_tuple(self, width, height, probability):
        self.whp_obj.add_tuple(width, height, probability)

    def generate_norm_distribution(self, h_mean, h_min, h_max, w_mean, w_min, w_max, std_h = 1, std_w = 1):
        a_w, b_w = (w_min - w_mean) / std_w, (w_max - w_min) / std_w
        a_h, b_h = (h_min - h_mean) / std_h, (h_max - h_mean) / std_h
        wrv = truncnorm.rvs(a_w, b_w, loc=w_mean, scale=std_w, size=10000)
        hrv = truncnorm.rvs(a_h, b_h, loc=h_mean, scale=std_h, size=10000)
        hist, xedges, yedges = np.histogram2d(wrv, hrv, bins=5, range=[[w_min, w_max], [h_min, h_max]])
        hist = hist / 10000
        w_vals = xedges[:-1] + 0.25
        h_vals = yedges[:-1] + 0.25
        for i, w_val in enumerate(w_vals):
            for j, h_val in enumerate(h_vals):
                self.add_tuple(w_val, h_val, hist[i,j])
    def generate_single_obstacle(self, w_val, h_val):
        self.add_tuple(w_val, h_val, 1)


class ShadowingRun:
    def __init__(self, x_min, x_max, y_min, y_max, height, whp, num_x=100, num_y=100):
        self.x = np.linspace(x_min, x_max, num_x)
        self.y = np.linspace(y_min, y_max, num_y)
        self.xy = np.zeros((self.x.shape[0], self.y.shape[0]))
        self.height = height
        self.whp = whp
    def __call__(self, x_t, y_t, z_t):
        for i, x in enumerate(self.x):
            for j, y in enumerate(self.y):
                logger.info('x: {}, y: {}'.format(x, y))
                positions = [{"xt": x_t, "yt": y_t, "xr": x, "yr": y, "zi": z_t, "zj": self.height}]
                whp_bind = self.whp.whp_obj.to_bind(0, max(y_t, y), 3, 5)

                sh_param = cm.Shadowing_Parameters(positions, whp_bind)
                sh_param_bind = sh_param.to_bind()
                self.xy[i,j] = chb.calculate_expectancy(sh_param_bind)
                logger.info('expected: {}'.format(self.xy[i,j]))
        return self.xy
    def get_exponential_form(self, intensity, time):
        return np.exp(-intensity*self.xy*time)
    def just_one_run(self, x_t, y_t, z_t, x, y, height = None):
        logger.info('x: {}, y: {}'.format(x, y))
        if height is None:
            height = self.height
        positions = [{"xt": x_t, "yt": y_t, "xr": x, "yr": y, "zi": z_t, "zj": height}]
        whp_bind = self.whp.whp_obj.to_bind(0.5, max(y_t, y), 3, 5)

        sh_param = cm.Shadowing_Parameters(positions, whp_bind)
        sh_param_bind = sh_param.to_bind()
        logger.info('expected: {}'.format(chb.calculate_expectancy(sh_param_bind)))


if __name__ == "__main__":
    wrap = WHP_Wrapper()
    h_min = 1
    w_min = 0.5
    h_max = 3
    w_max = 3
    h_mean = 1.5
    w_mean = 2
    std_h = 0.5
    std_w = 0.5
    w_single = 1.5
    h_single = 1.8
    wrap.generate_norm_distribution(h_mean, h_min, h_max, w_mean, w_min, w_max, std_h, std_w)
    #wrap.generate_single_obstacle(w_single,h_single)
    x_max = 5
    x_min = 0
    y_max = 3
    y_min = 0.51
    height = 1.8
    shadowing_instance = ShadowingRun(x_min, x_max, y_min, y_max, height, wrap)
    xled = 2
    #shadowing_instance.just_one_run(xled, 0.5, 3, 5, 3, 1)
    #shadowing_instance.just_one_run(xled, 0.5, 3, 5, 3, 1.6)
    #shadowing_instance.just_one_run(xled, 0.5, 3, 5, 3, 1.8)
    xy = shadowing_instance(xled,0.5,3)
    #print(xy)
    xy.tofile('at_{}_user_{}_normd.csv'.format(xled, height), sep=",")
