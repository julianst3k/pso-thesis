import ChannelSimulation.RunModel as rm
import numpy as np
import argparse
import simulation_runner as sr
from get_statistics import from_response_to_drms_and_power, average_by_location, statistics_by_transmitter, power_cummulative_distribution
from plot_statistics import plot_power_dBm
from yaml import load, dump, Loader
parser = argparse.ArgumentParser("")
parser.add_argument('--response', dest='response', type=str, help="The type of response that you want")
args = parser.parse_args()
config = 'config.yaml'
beta_arr = [20,25,35,40,50,55,65,70]
for beta in beta_arr:
    with open("config.yaml", "r+") as f:
        doc = load(f, Loader = Loader)
        doc["beta_r"] = beta
    with open("config.yaml", "w") as f:
        dump(doc, f)
    run_channel = rm.ModelRunWrapper(config, rm.Runner.CHANNEL)
    X = np.linspace(0.5,4.5,50)
    Y = np.linspace(0.5,2.5,25)
    #plot_power_dBm(X, Y, 0)
    sr.run_responses(X, Y, config)
    from_response_to_drms_and_power(X, Y, beta)
    average_by_location(X, Y, beta)
    power_cummulative_distribution(X, Y, beta)
