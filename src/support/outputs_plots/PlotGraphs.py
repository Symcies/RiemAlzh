import numpy

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import os
import time
import math

def PlotFinalOutputUnivariate():
    f = open("LastLineFile", 'r')
    last_line = f.readline().split()
    t0 = float(last_line[3]) #TauMean
    p  = float(last_line[2]) #P
    v0 = math.exp(float(last_line[5])) #exp(KsiMean)

    X = numpy.linspace(50,100,50)
    Y = []
    for x in X:
        Y.append(univarF(x, p, v0, t0))
    plt.figure()
    plt.plot(X, Y)

    plt.show()

def univarF(t,p,v0,t0):
    return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))

def PlotMultivarOutputWhileComputing(filename):
    # Activates interactivity
    plt.ion()

    # Shapes the figure: We have 2 rows, 3 columns
    fig = plt.figure()
    subplotNoise = fig.add_subplot(231)
    subplotNoise.set_title("Noise")
    subplotG = fig.add_subplot(234)
    subplotG.set_title("G")

    subplotMeanTau = fig.add_subplot(232)
    subplotMeanTau.set_title("Mean Tau")

    subplotVarTau = fig.add_subplot(235)
    subplotVarTau.set_title("Var Tau")

    subplotMeanKsi = fig.add_subplot(233)
    subplotMeanKsi.set_title("Mean Ksi")

    subplotVarKsi = fig.add_subplot(236)
    subplotVarKsi.set_title("Var Ksi")
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)

    line_vec = []
    while os.path.isfile(filename):
        #plt.clf()
        x_vec        = []
        noise        = []
        g_vec        = []
        tau_var_vec  = []
        tau_mean_vec = []
        ksi_var_vec  = []
        ksi_mean_vec = []
        for line in open(filename, 'r'):
            line_vec = line.split() #To change if lots of params?
            x_vec.append(line_vec[0])
            noise.append(line_vec[1])
            g_vec.append(line_vec[2])
            tau_mean_vec.append(line_vec[3])
            tau_var_vec.append(line_vec[4])
            ksi_mean_vec.append(line_vec[5])
            ksi_var_vec.append(line_vec[6])

        subplotNoise.plot(x_vec, noise, 'r')
        plt.pause(0.0001)
        subplotG.plot(x_vec, g_vec, 'r')
        plt.pause(0.0001)
        subplotMeanTau.plot(x_vec, tau_mean_vec, 'b')
        plt.pause(0.0001)
        subplotVarTau.plot(x_vec, tau_var_vec, 'b')
        plt.pause(0.0001)
        subplotMeanKsi.plot(x_vec, ksi_mean_vec, 'g')
        plt.pause(0.0001)
        subplotVarKsi.plot(x_vec, ksi_var_vec, 'g')
        plt.pause(0.0001)

    plt.ioff()
    plt.show(block = True)

    f = open("LastLineFile", 'w+')
    f.write(line_vec)

def PlotUnivarOutputWhileComputing(filename):
    # Activates interactivity
    plt.ion()

    # Shapes the figure: We have 2 rows, 3 columns
    fig = plt.figure()
    subplotNoise = fig.add_subplot(231)
    subplotNoise.set_title("Noise")
    subplotP = fig.add_subplot(234)
    subplotP.set_title("P")

    subplotMeanTau = fig.add_subplot(232)
    subplotMeanTau.set_title("Mean Tau")

    subplotVarTau = fig.add_subplot(235)
    subplotVarTau.set_title("Var Tau")

    subplotMeanKsi = fig.add_subplot(233)
    subplotMeanKsi.set_title("Mean Ksi")

    subplotVarKsi = fig.add_subplot(236)
    subplotVarKsi.set_title("Var Ksi")
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)

    line = ""
    while os.path.isfile(filename):
        #plt.clf()
        x_vec        = []
        noise        = []
        p_vec        = []
        tau_var_vec  = []
        tau_mean_vec = []
        ksi_var_vec  = []
        ksi_mean_vec = []
        for line in open(filename, 'r'):
            line_vec = line.split() #To change if lots of params?
            x_vec.append(line_vec[0])
            noise.append(line_vec[1])
            p_vec.append(line_vec[2])
            tau_mean_vec.append(line_vec[3])
            tau_var_vec.append(line_vec[4])
            ksi_mean_vec.append(line_vec[5])
            ksi_var_vec.append(line_vec[6])

        subplotNoise.plot(x_vec, noise, 'r')
        plt.pause(0.0001)
        subplotP.plot(x_vec, p_vec, 'r')
        plt.pause(0.0001)
        subplotMeanTau.plot(x_vec, tau_mean_vec, 'b')
        plt.pause(0.0001)
        subplotVarTau.plot(x_vec, tau_var_vec, 'b')
        plt.pause(0.0001)
        subplotMeanKsi.plot(x_vec, ksi_mean_vec, 'g')
        plt.pause(0.0001)
        subplotVarKsi.plot(x_vec, ksi_var_vec, 'g')
        plt.pause(0.0001)

        time.sleep(1)

    plt.ioff()
    plt.show(block = True)

    f = open("LastLineFile", 'w')
    f.write(line)
    f.close()

    return