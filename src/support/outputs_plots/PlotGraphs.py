import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import os
import time
import math

def PlotFinalOutputUnivariate(filename):
    # Extraction of variables
    f = open(filename, 'r')
    labels_aver = f.readline().split() #Unused at the moment
    aver_line = f.readline().split()
    labels_var = f.readline().split() #Unused at the moment
    tau_line = []
    ksi_line = []
    for line in f:
        l = line.split()
        tau_line.append(l[0])
        ksi_line.append(l[1])

    # Creation of X values
    X = numpy.linspace(40,110,70)
    plt.figure()

    # Average variables
    p  = float(aver_line[2]) #P
    aver_t0 = float(aver_line[3]) #TauMean
    aver_v0 = math.exp(float(aver_line[5])) #exp(KsiMean)

    # Mean
    aver_Y = []
    for x in X:
        aver_Y.append(getF("Univariate", x, p, aver_v0, aver_t0))
    plt.plot(X, aver_Y, 'r', linewidth=1)

    # Individual reals
    for i in range(len(tau_line)):
        Y = []
        t0 = float(tau_line[i]) #TauMean
        v0 = math.exp(float(ksi_line[i])) #exp(KsiMean)
        for x in X:
            Y.append(getF("Univariate", x, p, v0, t0))
        plt.plot(X, Y, linewidth=0.1)

    plt.show()


def PlotAndSelectPatientCurvesWithData(filename, map_list, type):
    print "in plot"
    values = PlotPatientCurves(filename, map_list, type, False)

    def update(patient_id):
        values[0][patient_id].set_visible(not values[0][patient_id].get_visible())
        plt.pause(0.0001)
        values[1][patient_id].set_visible(not values[1][patient_id].get_visible())
        plt.pause(0.0001)

    #plt.draw()

    while True:
        value = input("Get patient?")
        try:
            int(str(value), 10) + 1 #TODO: fix problem with octal numbers + add error management
        except Error:
            print "This is not a number, exiting the program."
            return
        update(int(value))

    plt.ioff()


def PlotPatientCurves(filename, map_list, type, visible):
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
        aver_Y.append(getF(type, x, p, aver_v0, aver_t0))
    splot.plot(X, aver_Y, 'r', linewidth=1, visible = True)

    list_lines = []
    list_points = []
    for i in range(len(tau_line)):
        # Individual reals
        Y = []
        t0 = float(tau_line[i]) #TauMean
        v0 = math.exp(float(ksi_line[i])) #exp(KsiMean)
        for x in X:
            Y.append(getF(type, x, p, v0, t0))
        line, = splot.plot(X, Y, linewidth=0.5, visible = visible)

        map = map_list[i]
        real_X = map.keys()
        real_Y = map.values()
        points, = splot.plot(real_X, real_Y, color = line.get_color(), linestyle = ' ', marker = '+', visible = visible)
        list_lines.append(line)
        list_points.append(points)

    return [list_lines, list_points]

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

def getF(type, t, p, v0, t0):
    if type == "Univariate":
        return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))
    if type == "Multivariate":
        return
    print type