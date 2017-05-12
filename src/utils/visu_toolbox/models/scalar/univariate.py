import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import utils.input_management as im


matplotlib.use('tkagg')

fig = plt.figure()
splot = fig.add_subplot(111)


def plot_curves(filename_pop, filename_ind, observations):
    plt.ion()

    pop_labels = im.read_population_parameters(filename_pop)
    number_param, param_dict = im.read_individual_parameters(filename_ind)

    plot_mean_curve(pop_labels)
    return plot_indiv_curves(param_dict, observations, pop_labels)


def f(t, p, v0, t0):
    return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))


def plot_mean_curve(pop_labels):
    # Creation of X values
    X = numpy.linspace(40, 110, 70)
    # Mean
    aver_Y = []
    for x in X:
        aver_Y.append(f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0]), pop_labels["taumean"][0]))
    splot.plot(X, aver_Y, 'r', linewidth=1, visible = True)


def plot_indiv_curves(param_dict, observations, pop_labels):
    color = numpy.random.rand(1, 3)

    # Individual curves
    X = numpy.linspace(40, 110, 70)
    list_lines = []
    list_points = []

    for i in range(len(param_dict["id"][0])): #For patient i
        # Individual reals
        t0 = param_dict["tau"][0][i] #TauMean
        v0 = math.exp(param_dict["ksi"][0][i]) #exp(KsiMean)

        Y = []

        for x in X:
            Y.append(f(x, pop_labels["p"][0], v0, t0))
        line, = splot.plot(X, Y, color = color[0], linewidth=0.5, visible = False)
        list_lines.append(line)

        real_X = observations[param_dict["id"][0][i]].obs.keys()
        real_Y = []
        for val in observations[param_dict["id"][0][i]].obs.values():
            real_Y.append(val[0])
        point, = splot.plot(real_X, real_Y, color = color[0], linestyle = ' ', marker = '+', visible = False)
        list_points.append(point)

    return list_lines, list_points