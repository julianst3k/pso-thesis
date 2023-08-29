import pandas as pd
import os
def read_pandas(pd_route, position):
    input_frame = pd.read_csv(pd_route)
    new_df = []
    for index, row in input_frame.iterrows():
        value = row["value"]
        particle = row["particle"]
        iteration = row["iter"]
        dictionary = {}
        dictionary["iter"] = iteration
        dictionary["sinr_h1"], dictionary["sinr_h2"] = value
        dictionary["betat"], dictionary["betar"], dictionary["fov"], dictionary["alphai"], dictionary["alphao"] = particle
        dictionary["x"], dictionary["y"] = position
        new_df.append(dictionary)
    new_pandas = pd.DataFrame(new_df)
    return new_pandas

def concat_pandas(pds):
    return pd.concat(pds)

if __name__ == "__main__":
    dataframe_aggregator = []
    for file in os.scandir("."):
        if file.endswith(".csv"):
            split_file = file.split("_")
            direction = split_file[-1]
            x = float(direction[0:2])
            y = float(direction[2:4])
            position_dataframe = pd.read_csv(file)
            dataframe_aggregator.append(position_dataframe)
    output_dataframe = concat_pandas(dataframe_aggregator)
    output_dataframe.to_csv("pso_opt.csv")
