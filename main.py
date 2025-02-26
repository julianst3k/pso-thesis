import ChannelSimulation.RunModel as rm
import numpy as np
import argparse
import simulation_runner as sr
from get_statistics import from_response_to_drms_and_power, average_by_location, statistics_by_transmitter
from plot_statistics import plot_power_dBm
parser = argparse.ArgumentParser("")
parser.add_argument('--response', dest='response', type=str, help="The type of response that you want")
args = parser.parse_args()
config = 'config.yaml'
run_channel = rm.ModelRunWrapper(config, rm.Runner.CHANNEL)

X = np.linspace(0.5,4.5,50)
Y = np.linspace(0.5,2.5,25)
#plot_power_dBm(X, Y, 0)
sr.run_responses(X, Y, config)