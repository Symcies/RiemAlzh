import numpy
import math
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('tkagg')

# brown1, dark green, light blue, purple,
# dark blue,  pink, black, grey,  orange, light green
color_map = ['#873D00' , '#007515', '#2436FF',
             '#00AFFA', '#9900FF', '#DB57FF', '#000000', '#808080',  '#FF7300', '#00D426']

class Multivariate:
    X = numpy.linspace(40, 110, 70)
    list_lines = []
    list_points = []

    def __init__(self, pop_params, indiv_param, number_param, observations, splot):
        self.aver_lines = []
        self.std_dev_lines = {"tau": [], "ksi": []}
        self.init_mean_curve(pop_params, number_param, splot)
        self.init_indiv_curves(indiv_param, observations, pop_params, number_param, splot)

    def init_mean_curve(self, pop_labels, number_param, splot):
        # Mean
        aver_Y = [[] for _ in xrange(int(number_param["w"]))]
        aver_tau = [[[],[]]  for _ in xrange(int(number_param["w"]))]
        aver_ksi = [[[],[]]  for _ in xrange(int(number_param["w"]))]
        for k in range(int(number_param["w"])):
            for x in self.X:
                aver_Y[k].append(
                    self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0, math.exp(pop_labels["ksimean"][0]),
                      pop_labels["taumean"][0]))
                aver_tau[k][0].append(self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0,math.exp(pop_labels["ksimean"][0]),
                                          pop_labels["taumean"][0] - math.sqrt(pop_labels["tauvar"][0])))
                aver_tau[k][1].append(self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0, math.exp(pop_labels["ksimean"][0]),
                                          pop_labels["taumean"][0] + math.sqrt(pop_labels["tauvar"][0])))
                aver_ksi[k][0].append(self.f(x, pop_labels["g"][0], pop_labels["deltas"][k], 0,
                                          math.exp(pop_labels["ksimean"][0] - math.sqrt(pop_labels["ksivar"][0])),
                                          pop_labels["taumean"][0]))
                aver_ksi[k][1].append(self.f(x, pop_labels["g"][0] , pop_labels["deltas"][k], 0,
                                          math.exp(pop_labels["ksimean"][0] + math.sqrt(pop_labels["ksivar"][0])),
                                          pop_labels["taumean"][0]))

            self.aver_lines.append(splot.plot(self.X, aver_Y[k], color=color_map[k], linewidth=1, visible=True))
            self.std_dev_lines["tau"].append([
                splot.plot(self.X, aver_tau[k][0], 'g', linewidth=1, linestyle='--', visible=False),
                splot.plot(self.X, aver_tau[k][1], 'g', linewidth=1, linestyle='--', visible=False)])
            self.std_dev_lines["ksi"].append([
                splot.plot(self.X, aver_ksi[k][0], 'g', linewidth=1, linestyle=':', visible=False),
                splot.plot(self.X, aver_ksi[k][1], 'g', linewidth=1, linestyle=':', visible=False)])

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
                line, = splot.plot(self.X, Y, color = color_map[k], linewidth=1, visible = False)
                lines.append(line)

                real_Y_k = []
                for val in observations[param_dict["id"][0][i]].obs.values():
                    real_Y_k.append(val[k])
                point, = splot.plot(real_X, real_Y_k, color = color_map[k], linestyle = ' ', marker = '+', visible = False)

                points.append(point)

            self.list_lines.append(lines)
            self.list_points.append(points)

    def f(self, t, g, delta, w, v0, t0):
        return 1 / (1 + g*math.exp((-w*math.pow(g*math.exp(-delta)+1,2)/g*math.exp(-delta)) - delta - v0*(t - t0)))

    def plot_mean(self):
        for i in range(len(self.aver_lines)):
            self.aver_lines[i][0].set_visible(not self.aver_lines[i][0].get_visible())
            plt.pause(0.0001)


    def plot_stand_dev(self, param):
        for i in range(len(self.std_dev_lines[param])):
            self.std_dev_lines[param][i][0][0].set_visible(not self.std_dev_lines[param][i][0][0].get_visible())
            self.std_dev_lines[param][i][1][0].set_visible(not self.std_dev_lines[param][i][1][0].get_visible())
            print self.std_dev_lines[param][i][1][0].get_visible()

    def plot_patients(self, list_patients_id, with_obs):
        for i in list_patients_id:
            for j in range(len(self.aver_lines)):
                self.list_lines[i][j].set_visible(not self.list_lines[i][j].get_visible())
                self.list_points[i][j].set_visible(with_obs)
