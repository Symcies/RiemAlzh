import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import math

def PlotPatientCurvesMultivariate(filename, isVisible, hasObservations, observationsMap):
    plt.ion()
    # Extraction of variables
    f = open(filename, 'r')
    labels_pop = f.readline().split() #Unused at the moment
    aver_line = f.readline().split()

    # Average variables
    aver = {}
    for i in range(len(labels_pop)):
        aver[labels_pop[i]] = float(aver_line[i])

    labels_ind = f.readline().split()
    # Average variables
    index_labels = {}
    for i in range(len(labels_ind)):
        index_labels[labels_ind[i]] = i
    id_line = []
    tau_line = []
    ksi_line = []
    w_lines = [[] for _ in xrange(int(aver["NumReal"]))]
    for line in f:
        l = line.split()
        id_line.append(float(l[index_labels["id"]]))
        tau_line.append(float(l[index_labels["Tau"]]))
        ksi_line.append(float(l[index_labels["Ksi"]]))
        for i in range(len(w_lines)):
            w_lines[i].append(float(l[index_labels["W"+i]]))

    # Creation of X values
    X = numpy.linspace(50,110,90)
    fig = plt.figure()
    splot = fig.add_subplot(111)

    # Mean
    aver_Y = [[] for _ in xrange(int(aver["NumReal"]))]
    colors = []
    for k in range(int(aver["NumReal"])):
        color = numpy.random.rand(1,3)
        colors.append(color)
        for x in X:
            aver_Y[k].append(fMultivariate(x, aver["G"], aver["Delta" + str(k)], 0, math.exp(aver["KsiMean"]), aver["TauMean"]))
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
            real_X = observationsMap[id_line[i]].keys()

        for k in range(int(aver["NumReal"])): #For curve k
            Y = []
            for x in X:
                Y.append(fMultivariate(x, aver["G"], aver["Delta" + str(k)], w_lines[k][i], v0, t0))
            line, = splot.plot(X, Y, color = colors[k][0], linewidth=0.5, visible = isVisible)
            lines.append(line)


            if hasObservations:
                real_Y_k = []
                for val in observationsMap[id_line[i]].values():
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
