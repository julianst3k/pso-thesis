import pandas as pd
import os
import re
import matplotlib.pyplot as plt
def read_pandas(pd_route, position):
    input_frame = pd.read_csv(pd_route)
    print(pd_route)
    new_df = []
    for index, row in input_frame.iterrows():
        value = row["value"]
        particle = row["particles"]
        iteration = row["iter"]
        dictionary = {}
        dictionary["iter"] = iteration
        string_value = re.sub("([\]|\[])", "", value)
        string_particle = re.sub("([\]|\[])", "", particle)
        vec_value = string_value.split(" ")
        particle_value = string_particle.split(" ")
        value = [float(string.replace(",", "")) for string in vec_value if string != '']
        particle = [float(string.replace(",", "")) for string in particle_value if  string != '']
        dictionary["sinr_h1"], dictionary["sinr_h2"] = value[0], value[1]
        dictionary["betat"], dictionary["betar"], dictionary["fov"], dictionary["alphao"], dictionary["alphai"] = particle[0], particle[1], particle[2], particle[3], particle[4]
        dictionary["x"], dictionary["y"] = position
        dictionary["index_value"] = min(value[0]/value[1], value[1]/value[0])
        new_df.append(dictionary)
    new_pandas = pd.DataFrame(new_df)
    return new_pandas

def concat_pandas(pds):
    return pd.concat(pds)
def process():
    dataframe_aggregator = []
    for file in os.scandir("."):
        if file.name.endswith(".csv"):
            split_file = file.name.split("_")
            direction = split_file[-1]
            try:
                x = float(direction[0:2])
                y = float(direction[2:4])
                position_dataframe = read_pandas(file.name, [x / 10, y / 10])
                dataframe_aggregator.append(position_dataframe)
            except ValueError:
                continue

    output_dataframe = concat_pandas(dataframe_aggregator)
    output_dataframe.to_csv("pso_opt.csv")
def filter_by_iter(file_name):
    unfiltered_dataframe = pd.read_csv(file_name)
    unfiltered_dataframe_sorted = unfiltered_dataframe.where(unfiltered_dataframe["iter"]==59).dropna().sort_values("index_value", ascending=False).groupby(["x","y"]).head(10)
    return unfiltered_dataframe_sorted
def filter_by_value():
    unfiltered_dataframe = pd.read_csv("pso_opt.csv")
    sinr2_filt = unfiltered_dataframe.sinr_h2 < -0.9
    sinr1_filt = unfiltered_dataframe.sinr_h1 < -0.9

    return unfiltered_dataframe.where(sinr1_filt & sinr2_filt).dropna()
def filter_by_iter_and_value(file_name):
    unfiltered_dataframe = pd.read_csv(file_name)
    iter_filt = unfiltered_dataframe["iter"]==59
    columns = ["betat","betar","fov","alphai","alphao"]
    quantiles = {}
    for column in columns:
        upp_quantile = unfiltered_dataframe[column].quantile(.75)
        low_quantile = unfiltered_dataframe[column].quantile(.25)
        quantiles[column] = {"low": low_quantile, "high": upp_quantile}
    for column in columns:
        upp_quantile_cond = unfiltered_dataframe[column] <= quantiles[column]["high"]
        low_quantile_cond = unfiltered_dataframe[column] >= quantiles[column]["low"]
        unfiltered_dataframe = unfiltered_dataframe[upp_quantile_cond & low_quantile_cond].dropna()
    unfiltered_dataframe_sorted = unfiltered_dataframe.where(unfiltered_dataframe["iter"]==59).dropna().sort_values("index_value", ascending=False).groupby(["x","y"]).head(10)
    return unfiltered_dataframe_sorted
def lambda_app():
    unfiltered_dataframe = pd.read_csv("pso_opt_iterval.csv")
    value_frame = unfiltered_dataframe[["sinr_h1", "sinr_h2"]]
    min_frame = value_frame.min(axis=1)
    max_frame = value_frame.max(axis=1)
    value_frame["min"] = min_frame
    value_frame["max"] = max_frame
    return_frame = value_frame.assign(Coef = lambda x: 1/((abs(x["sinr_h1"])+abs(x["sinr_h1"]))*x["min"]/x["max"]))
    frame_sum = return_frame["Coef"].sum()
    return_frame = return_frame.assign(Weight = lambda x: x["Coef"]/frame_sum)
    unfiltered_dataframe["Weights"] = return_frame["Weight"]
    column = ["betat", "betar", "fov", "alphai", "alphao"]
    for col in column:
        weight_col = unfiltered_dataframe[col]*unfiltered_dataframe["Weights"]
        unfiltered_dataframe["Weighted {}".format(col)] = weight_col
    return unfiltered_dataframe
def sum_values():
    column = ["betat", "betar", "fov", "alphai", "alphao"]
    unfiltered_dataframe = pd.read_csv("values.csv")

    for col in column:
        weight_string = "Weighted {}".format(col)
        sum_weight = unfiltered_dataframe[weight_string].sum()
        print("{col} Weighted Average is {avg}".format(col=col, avg=sum_weight))
def process_three(file_name):
    read_frame = pd.read_csv(file_name)
    hist = read_frame.hist(column=["betat","betar","fov","alphao","alphai"], bins=25, layout=(1,5), grid=False, color='maroon', rwidth=0.9)
    print(hist)
    title = ["(a)","(b)","(c)","(d)","(e)"]
    x_label = "Angle [°]"
    y_label = "Frequency"
    for i, x_np in enumerate(hist):
        for j, ax in enumerate(x_np):
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off",
                          labelleft="on")
            if j==0:
                ax.set_ylabel(y_label)
            ax.set_xlabel(x_label)
            ax.set_title(title[j])
            ax.plot()
    print(read_frame[["betat","betar","fov","alphai","alphao"]].describe())
    plt.show()
def plot_evolution(file_name):
    read_frame = pd.read_csv(file_name)
    interest_variable = [.3,.5]
    xfilter = read_frame["x"] == interest_variable[0]
    yfilter = read_frame["y"] == interest_variable[1]
    filtered_frame = read_frame.where(xfilter & yfilter).dropna()
    filtered_frame['sinr_h1'] = filtered_frame['sinr_h1'].apply(lambda x: -x)
    filtered_frame['sinr_h2'] = filtered_frame['sinr_h2'].apply(lambda x: -x)
    print(filtered_frame)
    first_iteration = filtered_frame.where(filtered_frame["iter"]==0)
    tenth_iteration = filtered_frame.where(filtered_frame["iter"]==20)
    twenty_iteration = filtered_frame.where(filtered_frame["iter"]==40)
    sixty_iteration = filtered_frame.where(filtered_frame["iter"]==59)
    iteration_agg = [first_iteration, tenth_iteration, twenty_iteration, sixty_iteration]
    colores = ['DarkBlue', 'Teal', 'Maroon', 'Black']
    iteration_order = [0, 20, 40, 60]
    for i, iteration in enumerate(iteration_agg):
        values = iteration[["sinr_h1", "sinr_h2"]]
        if i == 0:
            scatter = values.plot.scatter(x='sinr_h1', y = 'sinr_h2', c = colores[i])
            values.plot.line(x='sinr_h1', y='sinr_h2', c=colores[i], ax = scatter, label = "Iteration {}".format(iteration_order[i]))
        else:
            values.plot.scatter(x='sinr_h1', y = 'sinr_h2', c = colores[i], ax = scatter)
            values.plot.line(x='sinr_h1', y='sinr_h2', c=colores[i], ax = scatter, label = "Iteration {}".format(iteration_order[i]))

        scatter.set_xlabel("SINR (Symbol 1)")
        scatter.set_ylabel("SINR (Symbol 2)")
    plt.show()



    title = ["(a)","(b)","(c)","(d)","(e)"]
    x_label = "Angle [°]"
    y_label = "Frequency"
    for i, x_np in enumerate(hist):
        for j, ax in enumerate(x_np):
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off",
                          labelleft="on")
            if j==0:
                ax.set_ylabel(y_label)
            ax.set_xlabel(x_label)
            ax.set_title(title[j])

            ax.plot()
    print(read_frame[["betat","betar","fov","alphai","alphao"]].describe())
    plt.show()

if __name__ == "__main__":
    """
    dir = os.scandir()
    dataframes = []
    for file in dir:
        filename = file.name
        if filename.endswith(".csv") and "_final_" in filename:
            position = filename.split("_")[-1]
            try:
                coord_x, coord_y = float(position[0:2])/10, float(position[2:4])/10
                file_pandas = read_pandas(filename, [coord_x, coord_y])
                dataframes.append(file_pandas)
            except ValueError:
                continue
    final_dataframe = pd.concat(dataframes)
    """
    process_three("filtered_pso.csv")
