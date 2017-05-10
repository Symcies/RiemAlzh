import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

import os
import time

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