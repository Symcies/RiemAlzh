import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt


# brown1, dark green, light blue, purple,
# dark blue,  pink, black, grey,  orange, light green
color_map = ['#873D00' , '#007515', '#2436FF',
             '#00AFFA', '#9900FF', '#DB57FF', '#000000', '#808080',  '#FF7300', '#00D426']

import models.scalar.multivariate as multivariate
import models.scalar.univariate as univariate
import utils.input_management as obs
import utils.math


def plot(output_pop_path, output_ind_path, observations_path, type):
    fig = plt.figure()
    splot = fig.add_subplot(111)

    list_lines, list_points = plot_whole_model(output_pop_path, output_ind_path, type, obs.extract_observations(observations_path), splot)
    plt.pause(3)

    while True:
        value = raw_input("Get patient? ")
        if utils.math.IsInt(value):
            if -1 < int(value) < len(list_lines):
                plot_patient(int(value), list_lines, list_points, type)
            else:
                print "The number of patient are in [0," + str(len(list_lines)-1) + "], and you entered " + value
        else:
            code = user_actions(str(value))
            if code == "all":
                change_curves_visibility(True, list_lines, list_points, type)
            elif code == "none":
                change_curves_visibility(False, list_lines, list_points, type)

    plt.plot(block = True)


def plot_whole_model(output_pop_path, output_ind_path, type, map_list, splot):
    if type == "Univariate":
        return univariate.plot_curves(output_pop_path, output_ind_path, map_list, splot)
    elif type == "Multivariate":
        return multivariate.plot_curves(output_pop_path, output_ind_path, map_list, splot)


def plot_patient(patient_id, list_lines, list_points, type):
    plt.ion()
    if type == "Multivariate":
        for i in range(4):
            list_lines[patient_id][i].set_visible(not list_lines[patient_id][i].get_visible())
            list_points[patient_id][i].set_visible(not list_points[patient_id][i].get_visible())
            plt.pause(0.0001)
    elif type == "Univariate":
        list_lines[patient_id].set_visible(not list_lines[patient_id].get_visible())
        list_points[patient_id].set_visible(not list_points[patient_id].get_visible())
        plt.pause(0.0001)

    plt.pause(0.0001)

    plt.ioff()


def change_curves_visibility(boolean, list_lines, list_points, type):
    plt.ion()
    if type == "Multivariate":
        for i in range(len(list_lines)):
            for j in range(4):
                list_lines[i][j].set_visible(boolean)
                list_points[i][j].set_visible(False)
    elif type == "Univariate":
        for i in range(len(list_lines)):
            list_lines[i].set_visible(boolean)
            list_points[i].set_visible(False)

    plt.pause(0.0001)
    plt.ioff()


def user_actions(value):
    if value == "hold" or value == "h" or value == "pause" or value == "p":
        duration = input("How long do you want to hold the display (in sec)?")
        if utils.math.IsInt(duration):
            plt.pause(int(duration))
        return 0
    if value == "exit" or value == "quit" or value == "q":
        exit(0)
    if value == "all" or value == "none":
        return value