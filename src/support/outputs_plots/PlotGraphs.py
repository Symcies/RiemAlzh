import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import os
import time
import math

def PlotFinalOutput(filename, type):
    lines = PlotPatientCurves(filename, [], type, False, False)

    plt.show(block = True)
    # for line in lines:
    #     line.set_visible(True)
    #     line.set_width(0.3)
    #     plt.pause(0.0001)

def PlotAndSelectPatientCurvesWithData(filename, map_list, type):
    values = PlotPatientCurves(filename, map_list, type, False, True)

    def update(patient_id):
        if type == "Multivariate":
            for i in range(4):
                values[0][patient_id][i].set_visible(not values[0][patient_id][i].get_visible())
                plt.pause(0.0001)
        elif type == "Univariate":
            values[0][patient_id].set_visible(not values[0][patient_id].get_visible())
            plt.pause(0.0001)

        values[1][patient_id].set_visible(not values[1][patient_id].get_visible())
        plt.pause(0.0001)

    while True:
        value = input("Get patient?")
        try:
            int(str(value), 10) + 1 #TODO: fix problem with octal numbers + add error management
        except Error:
            print "This is not a number, exiting the program."
            return
        update(int(value))

    plt.ioff()


def PlotOutputWhileComputing(filename):
    # Activates interactivity
    plt.ion()

    f = open(filename, 'r')
    labels = f.readline().split()
    f.close()

    fig = plt.figure()
    subplots = []
    for i in range(1, len(labels)):
        subplot =  fig.add_subplot(2,int(len(labels)/2),i)
        subplot.set_title(labels[i])
        subplots.append(subplot)

    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)

    while os.path.isfile(filename):
        f = open(filename, 'r')
        labels = f.readline().split()
        values = [[] for _ in xrange(len(labels))]
        for line in f:
            line_vec = line.split() #To change if lots of params?
            for i in range(len(labels)):
                values[i].append(line_vec[i])

        for i in range(1, len(labels)):
            subplots[i-1].cla()
            subplots[i-1].plot(values[0], values[i], 'b')
            subplots[i-1].set_title(labels[i])
            plt.pause(0.0001)

        time.sleep(1)

    plt.ioff()
    plt.show(block = True)

    return

def PlotPatientCurves(filename, map_list, type, isVisible, hasObservations):
    if type == "Univariate":
        return PlotPatientCurvesUnivariate(filename, map_list, isVisible, hasObservations)
    elif type == "Multivariate":
        return PlotPatientCurvesMultivariate(filename, map_list, isVisible, hasObservations)

def PlotPatientCurvesUnivariate(filename, observationsMap, isVisible, hasObservations):
    plt.ion()
    # Extraction of variables
    f = open(filename, 'r')
    f.readline().split() #Unused at the moment
    aver_line = f.readline().split()
    f.readline().split() #Unused at the moment
    tau_line = []
    ksi_line = []
    for line in f:
        l = line.split()
        tau_line.append(l[0])
        ksi_line.append(l[1])

    # Creation of X values
    X = numpy.linspace(40,110,70)
    fig = plt.figure()
    splot = fig.add_subplot(111)

    # Average variables
    p  = float(aver_line[2]) #P
    aver_t0 = float(aver_line[3]) #TauMean
    aver_v0 = math.exp(float(aver_line[5])) #exp(KsiMean)

    # Mean
    aver_Y = []
    for x in X:
        aver_Y.append(fUnivariate(x, p, aver_v0, aver_t0))
    splot.plot(X, aver_Y, 'r', linewidth=1, visible = True)

    list_lines = []
    list_points = []
    for i in range(len(tau_line)):
        # Individual reals
        Y = []
        t0 = float(tau_line[i]) #TauMean
        v0 = math.exp(float(ksi_line[i])) #exp(KsiMean)
        for x in X:
            Y.append(fUnivariate(x, p, v0, t0))
        line, = splot.plot(X, Y, linewidth=0.5, visible = isVisible)
        list_lines.append(line)

        if hasObservations:
            map = observationsMap[i]
            real_X = map.keys()
            real_Y = map.values()
            points, = splot.plot(real_X, real_Y, color = line.get_color(), linestyle = ' ', marker = '+', visible = isVisible)
            list_points.append(points)

    if hasObservations:
        return [list_lines, list_points]
    return list_lines

def PlotPatientCurvesMultivariate(filename, observationsMap, isVisible, hasObservations):
    plt.ion()
    # Extraction of variables
    f = open(filename, 'r')
    f.readline().split() #Unused at the moment
    aver_line = f.readline().split()
    f.readline().split() #Unused at the moment
    tau_line = []
    ksi_line = []
    w_lines = [[] for _ in xrange(int(aver_line[5]))]
    for line in f:
        l = line.split()
        tau_line.append(float(l[0]))
        ksi_line.append(float(l[1]))
        for i in range(len(w_lines)):
            w_lines[i].append(float(l[2+i]))

    # Creation of X values
    X = numpy.linspace(60,90,60)
    fig = plt.figure()
    splot = fig.add_subplot(111)

    # Average variables
    g  = float(aver_line[2]) #P
    aver_t0 = float(aver_line[3]) #TauMean
    aver_v0 = math.exp(float(aver_line[4])) #exp(KsiMean)
    deltas = []
    for i in range(int(aver_line[5])):
        deltas.append(float(aver_line[6 + i]))
    print deltas

    # Mean
    aver_Y = [[] for _ in xrange(int(aver_line[5]))]
    for k in range(len(aver_Y)):
        for x in X:
            aver_Y[k].append(fMultivariate(x, g, deltas[k], 0, aver_v0, aver_t0))
        splot.plot(X, aver_Y[k], 'r', linewidth=1, visible = True)

    list_lines = []
    list_points = []
    for i in range(len(tau_line)):
        # Individual reals
        Y = [[] for _ in xrange(int(aver_line[5]))]
        t0 = tau_line[i] #TauMean
        v0 = math.exp(ksi_line[i]) #exp(KsiMean)
        lines = []
        for k in range(len(Y)):
            for x in X:
                Y[k].append(fMultivariate(x, g, deltas[k], w_lines[k][i], v0, t0))
            line, = splot.plot(X, Y[k], color = line.get_color(), linewidth=0.5, visible = isVisible)
            lines.append(line)
        list_lines.append(lines)

        if hasObservations:
            map = observationsMap[i]
            real_X = map.keys()
            real_Y = map.values()
            points, = splot.plot(real_X, real_Y, color = line.get_color(), linestyle = ' ', marker = '+', visible = isVisible)
            list_points.append(points)

    if hasObservations:
        return [list_lines, list_points]
    return list_lines



def fUnivariate(t, p, v0, t0):
    return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))

def fMultivariate(t, g, delta, w, v0, t0):
    return 1 / (1 + g*math.exp((-w*math.pow(g*math.exp(-delta)+1,2)/g*math.exp(-delta)) - delta - v0*(t - t0)))