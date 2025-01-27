import ChannelSimulation.RunModel as rm
import numpy as np
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('--response', dest='response', type=str, help="The type of response that you want")
args = parser.parse_args()
config = 'config.yaml'
rm
run_response = rm.ModelRunWrapper(config, rm.Runner.RESPONSE)
run_channel = rm.ModelRunWrapper(config, rm.Runner.CHANNEL)

run_response.update([1,1,1.8])
run_channel.update([1,1,1.8])


response = np.array(run_response(args.response))
responselos = np.array(run_response("los"))
responsescat = np.array(run_response("scatt"))
responsenlos = np.array(run_response("nlos"))

response.tofile("response_file.csv", sep=",")
responselos.tofile("response_los_file.csv", sep=",")
responsescat.tofile("response_scat_file.csv", sep=",")
responsenlos.tofile("response_nlos_file.csv", sep=",")
