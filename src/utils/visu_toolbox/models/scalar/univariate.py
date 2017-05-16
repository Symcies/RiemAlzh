import numpy
import math
import matplotlib.pyplot as plt


class Univariate:
    X = numpy.linspace(40, 110, 70)
    list_lines = []
    list_points = []

    def __init__(self, pop_params, indiv_param, observations, splot):
        self.std_dev_lines = {}
        self.init_mean_curve(pop_params, splot)
        self.init_indiv_curves(indiv_param, observations, pop_params, splot)

    def init_mean_curve(self, pop_labels, splot):
        # Mean
        aver_Y = []
        aver_tau = [[],[]]
        aver_ksi = [[],[]]
        for x in self.X:
            aver_Y.append(self.f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0]), pop_labels["taumean"][0]))
            aver_tau[0].append(self.f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0]), pop_labels["taumean"][0] - math.sqrt(pop_labels["tauvar"][0])))
            aver_tau[1].append(self.f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0]), pop_labels["taumean"][0] + math.sqrt(pop_labels["tauvar"][0])))
            aver_ksi[0].append(self.f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0] - math.sqrt(pop_labels["ksivar"][0])), pop_labels["taumean"][0]))
            aver_ksi[1].append(self.f(x, pop_labels["p"][0], math.exp(pop_labels["ksimean"][0] + math.sqrt(pop_labels["ksivar"][0])), pop_labels["taumean"][0]))

        self.aver_line = splot.plot(self.X, aver_Y, 'r', linewidth=1, visible = False)
        self.std_dev_lines["tau"] = [splot.plot(self.X, aver_tau[0], 'g', linewidth=1, linestyle = '--', visible = False),
                                     splot.plot(self.X, aver_tau[1], 'g', linewidth=1, linestyle = '--', visible = False)]
        self.std_dev_lines["ksi"] = [splot.plot(self.X, aver_ksi[0], 'g', linewidth=1, linestyle = ':', visible = False),
                                     splot.plot(self.X, aver_ksi[1], 'g', linewidth=1, linestyle = ':', visible = False)]

    def init_indiv_curves(self, param_dict, observations, pop_labels, splot):
        color = numpy.random.rand(1, 3)

        # Individual curves

        for i in range(len(param_dict["id"][0])): #For patient i
            # Individual reals
            t0 = param_dict["tau"][0][i] #TauMean
            v0 = math.exp(param_dict["ksi"][0][i]) #exp(KsiMean)

            Y = []

            for x in self.X:
                Y.append(self.f(x, pop_labels["p"][0], v0, t0))
            line, = splot.plot(self.X, Y, color = color[0], linewidth=0.5, visible = False)
            self.list_lines.append(line)

            real_X = observations[param_dict["id"][0][i]].obs.keys()
            real_Y = []
            for val in observations[param_dict["id"][0][i]].obs.values():
                real_Y.append(val[0])
            point, = splot.plot(real_X, real_Y, color = color[0], linestyle = ' ', marker = '+', visible = False)
            self.list_points.append(point)

    def f(self, t, p, v0, t0):
        return 1/(1+(math.exp(-p))*math.exp(-v0*(t-t0)))

    def plot_mean(self):
        self.aver_line[0].set_visible(not self.aver_line[0].get_visible())

    def plot_stand_dev(self, param):
        self.std_dev_lines[param][0][0].set_visible(not self.std_dev_lines[param][0][0].get_visible())
        self.std_dev_lines[param][1][0].set_visible(not self.std_dev_lines[param][1][0].get_visible())

    def plot_patients(self, list_patients_id, with_obs):
        for i in list_patients_id:
            self.list_lines[i].set_visible(not self.list_lines[i].get_visible())
            self.list_points[i].set_visible(with_obs)