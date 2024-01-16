import ChannelSimulation.RunModel as rm
import numpy as np

config = 'config.yaml'

run_response = rm.ModelRunWrapper(config, rm.Runner.RESPONSE)
run_channel = rm.ModelRunWrapper(config, rm.Runner.CHANNEL)

run_response.update([1,1,1.8])
run_channel.update([1,1,1.8])

response = np.array(run_response())

print(response.shape)
response.tofile("response_file.csv", sep=",")