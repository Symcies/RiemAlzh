import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import math

def PlotPatientCurvesMultivariate(filename_pop, filename_ind, isVisible, hasObservations, observationsMap):
    plt.ion()
    pop_labels = read_population_parameters(filename_pop)

    # Extraction of average variables
    f_ind = open(filename_ind, 'r')
    param_nums = {}
    name_params = f_ind.readline().split()
    num_params = f_ind.readline().split()
    for i in range(len(name_params)):
        param_nums[name_params[i]] = num_params[i]

    # Individual variables dictionnary
    labels_ind = f_ind.readline().split()
    # Average variables
    index_labels = {}
    for i in range(len(labels_ind)):
        index_labels[labels_ind[i]] = i
    id_line = []
    tau_line = []
    ksi_line = []
    w_lines = [[] for _ in xrange(int(param_nums["W"]))]
    for line in f_ind:
        l = line.split()
        id_line.append(float(l[index_labels["id"]]))
        tau_line.append(float(l[index_labels["Tau"]]))
        ksi_line.append(float(l[index_labels["Ksi"]]))
        for i in range(len(w_lines)):
            w_lines[i].append(float(l[index_labels["W"+str(i)]]))

    # Creation of X values
    X = numpy.linspace(50,110,90)
    fig = plt.figure()
    splot = fig.add_subplot(111)

    # Mean
    aver_Y = [[] for _ in xrange(int(param_nums["W"]))]
    colors = []
    for k in range(int(param_nums["W"])):
        color = numpy.random.rand(1,3)
        colors.append(color)
        for x in X:
            aver_Y[k].append(fMultivariate(x, pop_labels["g"][0], pop_labels["deltas"][k], 0, math.exp(pop_labels["ksimean"][0]), pop_labels["taumean"][0]))
        splot.plot(X, aver_Y[k], color = colors[k][0], linewidth=1, visible = True)


    # Individual curves
    list_lines = []
    list_points = []

    for i in range(len(id_line)): #For patient i
        # Individual reals
        t0 = tau_line[i] #TauMean
        v0 = math.exp(ksi_line[i]) #exp(KsiMean)
        lines = []
        points = []

        if hasObservations:
            real_X = observationsMap[id_line[i]].obs.keys()

        for k in range(int(param_nums["W"])): #For curve k
            Y = []
            for x in X:
                Y.append(fMultivariate(x, pop_labels["g"][0], pop_labels["deltas"][k], w_lines[k][i], v0, t0))
            line, = splot.plot(X, Y, color = colors[k][0], linewidth=0.5, visible = isVisible)
            lines.append(line)


            if hasObservations:
                real_Y_k = []
                for val in observationsMap[id_line[i]].obs.values():
                    real_Y_k.append(val[k])
                point, = splot.plot(real_X, real_Y_k, color = colors[k][0], linestyle = ' ', marker = '+', visible = isVisible)
                points.append(point)

        list_lines.append(lines)
        list_points.append(points)


    if hasObservations:
        return [list_lines, list_points]
    return list_lines


def fMultivariate(t, g, delta, w, v0, t0):
    return 1 / (1 + g*math.exp((-w*math.pow(g*math.exp(-delta)+1,2)/g*math.exp(-delta)) - delta - v0*(t - t0)))


def read_population_parameters(file_name):
    params = {}
    f = open(file_name, 'r')

    for line in f:
        line = line.lower()

        # Extract parameters
        list_elements = line.split()
        param_name = list_elements[0]
        param_values = [float(list_elements[i]) for i in range(1, len(list_elements))]

        # Add to existing dict
        params[param_name] = param_values

    return params