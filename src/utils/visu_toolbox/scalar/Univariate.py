import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import math

def PlotPatientCurvesUnivariate(filename, isVisible, hasObservations, observationsMap):
    plt.ion()
    # Extraction of variables
    f = open(filename, 'r')
    labels = f.readline().split() #Unused at the moment
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
    aver = {}
    for i in range(len(labels)):
        aver[labels[i]] = float(aver_line[i])

    # Mean
    aver_Y = []
    for x in X:
        aver_Y.append(fUnivariate(x, aver["P"], math.exp(aver["KsiMean"]), aver["TauMean"]))
    splot.plot(X, aver_Y, 'r', linewidth=1, visible = True)

    list_lines = []
    list_points = []
    for i in range(len(tau_line)):
        # Individual reals
        Y = []
        t0 = float(tau_line[i]) #TauMean
        v0 = math.exp(float(ksi_line[i])) #exp(KsiMean)
        color = numpy.random.rand(1,3)
        for x in X:
            Y.append(fUnivariate(x, aver["P"], v0, t0))
        line, = splot.plot(X, Y, linewidth=0.5, visible = isVisible, color = color)
        list_lines.append(line)

        if hasObservations:
            map = observationsMap[i]
            real_X = map.keys()
            real_Y = []
            for val in map.values():
                real_Y.append(val[0])
            points, = splot.plot(real_X, real_Y, color = color, linestyle = ' ', marker = '+', visible = isVisible)
            list_points.append(points)

    if hasObservations:
        return [list_lines, list_points]
    return list_lines

def fUnivariate(t, p, v0, t0):
    return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))
