import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import utils.input_management as im

matplotlib.use('tkagg')


colors = []


def plot_curves(filename_pop, filename_ind, observations, splot):
    plt.ion()
    
    pop_labels = im.read_population_parameters(filename_pop)
    number_param, param_dict = im.read_individual_parameters(filename_ind)

    init_colors(int(number_param["w"]))

    plot_mean_curve(number_param, pop_labels, splot)
    return plot_indiv_curves(param_dict, observations, pop_labels, number_param, splot)


def f(t, g, delta, w, v0, t0):
    return 1 / (1 + g*math.exp((-w*math.pow(g*math.exp(-delta)+1,2)/g*math.exp(-delta)) - delta - v0*(t - t0)))


def plot_mean_curve(number_param, pop_labels, splot):
    # Creation of X values
    X = numpy.linspace(50, 110, 90)
    aver_Y = [[] for _ in xrange(int(number_param["w"]))]
    for k in range(int(number_param["w"])):
        for x in X:
            aver_Y[k].append(f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0, math.exp(pop_labels["ksimean"][0]),
                               pop_labels["taumean"][0]))

        splot.plot(X, aver_Y[k], color=colors[k][0], linewidth=1, visible=True)


def plot_indiv_curves(param_dict, observations, pop_labels, number_param, splot):
    # Individual curves
    X = numpy.linspace(50, 110, 90)
    list_lines = []
    list_points = []

    for i in range(len(param_dict["id"][0])): #For patient i
        # Individual reals
        t0 = param_dict["tau"][0][i] #TauMean
        v0 = math.exp(param_dict["ksi"][0][i]) #exp(KsiMean)
        lines = []
        points = []

        real_X = observations[param_dict["id"][0][i]].obs.keys()

        for k in range(int(number_param["w"])): #For curve k
            Y = []
            for x in X:
                Y.append(f(x, pop_labels["g"][0], pop_labels["deltas"][k], param_dict["w"][k][i], v0, t0))
            line, = splot.plot(X, Y, color = colors[k][0], linewidth=0.5, visible = False)
            lines.append(line)

            real_Y_k = []
            for val in observations[param_dict["id"][0][i]].obs.values():
                real_Y_k.append(val[k])
            point, = splot.plot(real_X, real_Y_k, color = colors[k][0], linestyle = ' ', marker = '+', visible = False)
            points.append(point)

        list_lines.append(lines)
        list_points.append(points)

    return list_lines, list_points


def init_colors(n):
    for k in range(n):
        color = numpy.random.rand(1, 3)
        colors.append(color)
