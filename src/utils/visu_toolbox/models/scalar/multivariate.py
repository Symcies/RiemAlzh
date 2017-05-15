import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import utils.input_management as im

matplotlib.use('tkagg')

colors = []


class Multivariate:
    X = numpy.linspace(40, 110, 70)
    list_lines = []
    list_points = []

    def __init__(self, pop_params, indiv_param, number_param, observations, splot):
        self.init_colors(int(number_param["w"]))
        self.aver_lines = []
        self.init_mean_curve(pop_params, number_param, splot)
        self.init_indiv_curves(indiv_param, observations, pop_params, number_param, splot)

    def init_mean_curve(self, pop_labels, number_param, splot):
        # Mean
        aver_Y = [[] for _ in xrange(int(number_param["w"]))]
        for k in range(int(number_param["w"])):
            for x in self.X:
                aver_Y[k].append(
                    self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0, math.exp(pop_labels["ksimean"][0]),
                      pop_labels["taumean"][0]))

            self.aver_lines.append(splot.plot(self.X, aver_Y[k], color=colors[k][0], linewidth=1, visible=True))

    def init_indiv_curves(self, param_dict, observations, pop_labels, number_param, splot):
        for i in range(len(param_dict["id"][0])): #For patient i
            # Individual reals
            t0 = param_dict["tau"][0][i] #TauMean
            v0 = math.exp(param_dict["ksi"][0][i]) #exp(KsiMean)
            lines = []
            points = []

            real_X = observations[param_dict["id"][0][i]].obs.keys()

            for k in range(int(number_param["w"])): #For curve k
                Y = []
                for x in self.X:
                    Y.append(self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], param_dict["w"][k][i], v0, t0))
                line, = splot.plot(self.X, Y, color = colors[k][0], linewidth=0.5, visible = False)
                lines.append(line)

                real_Y_k = []
                for val in observations[param_dict["id"][0][i]].obs.values():
                    real_Y_k.append(val[k])
                point, = splot.plot(real_X, real_Y_k, color = colors[k][0], linestyle = ' ', marker = '+', visible = False)
                points.append(point)

            self.list_lines.append(lines)
            self.list_points.append(points)

    def f(self, t, g, delta, w, v0, t0):
        return 1 / (1 + g*math.exp((-w*math.pow(g*math.exp(-delta)+1,2)/g*math.exp(-delta)) - delta - v0*(t - t0)))

    def plot_mean(self):
        for i in range(len(self.aver_lines)):
            self.aver_lines[i][0].set_visible(not self.aver_lines[i][0].get_visible())
            plt.pause(0.0001)

    def plot_patients(self, list_patients_id, with_obs):
        for i in list_patients_id:
            for j in range(len(self.aver_lines)):
                self.list_lines[i][j].set_visible(not self.list_lines[i][j].get_visible())
                self.list_points[i][0].set_visible(with_obs)

    def init_colors(self, n):
        for k in range(n):
            color = numpy.random.rand(1, 3)
            colors.append(color)
