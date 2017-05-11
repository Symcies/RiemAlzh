import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

from scalar import Univariate
from scalar import Multivariate
import utils

def PlotFinalOutput(output_pop_path, output_ind_path, observations_path, type):
    PlotAndSelectPatientCurvesWithData(output_pop_path, output_ind_path, type, utils.readInput(observations_path))

    plt.show(block = True)

def PlotAndSelectPatientCurvesWithData(output_pop_path, output_ind_path, type, map_list):
    values = PlotPatientCurves(output_pop_path, output_ind_path, type, False, True, map_list)
    plt.pause(3)

    def update(patient_id):
        plt.ion()
        if type == "Multivariate":
            for i in range(4):
                values[0][patient_id][i].set_visible(not values[0][patient_id][i].get_visible())
                values[1][patient_id][i].set_visible(not values[1][patient_id][i].get_visible())
                plt.pause(0.0001)
        elif type == "Univariate":
            values[0][patient_id].set_visible(not values[0][patient_id].get_visible())
            values[1][patient_id].set_visible(not values[1][patient_id].get_visible())
            plt.pause(0.0001)

        plt.pause(0.0001)
        plt.ioff()

    def updateAll(boolean):
        plt.ion()
        if type == "Multivariate":
            for i in range(len(values[0])):
                for j in range(4):
                    values[0][i][j].set_visible(boolean)
                    values[1][i][j].set_visible(False)
        elif type == "Univariate":
            for i in range(len(values[0])):
                values[0][i].set_visible(boolean)
                values[1][i].set_visible(False)

        plt.pause(0.0001)
        plt.ioff()

    while True:
        value = raw_input("Get patient?")
        if utils.IsInt(value):
            if int(value) < len(values[0]) and int(value)>-1:
                update(int(value))
            else:
                print "The number of patient are in [0," + str(len(values[0])-1) + "], and you entered " + value
        else:
            code = UserMessages(str(value))
            if code == 0: #quit
                break
            elif code == 2: #all
                updateAll(True)
            elif code == 3: #none
                updateAll(False)

    plt.plot(block = True)

def PlotPatientCurves(output_pop_path, output_ind_path, type, isVisible, hasObservations, map_list):
    if type == "Univariate":
        return Univariate.PlotPatientCurvesUnivariate(output_pop_path, output_ind_path, isVisible, hasObservations, map_list)
    elif type == "Multivariate":
        return Multivariate.PlotPatientCurvesMultivariate(output_pop_path, output_ind_path, isVisible, hasObservations, map_list)

def UserMessages(value):
    if value == "hold" or value == "h" or value == "pause" or value == "p":
        duration = input("How long do you want to hold the display (in sec)?")
        if utils.IsInt(duration):
            plt.pause(int(duration))
        return 1
    if value == "exit" or value == "quit" or value == "q":
        return 0
    if value == "all":
        return 2
    if value == "none":
        return 3
    return False