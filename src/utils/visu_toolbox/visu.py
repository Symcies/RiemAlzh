import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.io

matplotlib.use('tkagg')

from models.scalar.univariate import Univariate
from models.scalar.multivariate import Multivariate

from utils.input_management import read_individual_parameters, read_population_parameters, extract_observations
from models.networks.fast_network import *

class Visu:
    '''
    Documentation
    '''

    #model_type = ""

    def __init__(self, type, population_params, individual_params, observations):
        self.fig = plt.figure()
        self.splot = self.fig.add_subplot(111)

        self.model_type = type.lower()
        self.pop_params = read_population_parameters(population_params)
        self.number_param, self.indiv_params = read_individual_parameters(individual_params)
        self.observations = extract_observations(observations)
        self.model = self.init_model()


    def init_model(self):
        if self.model_type == "univariate":
            return Univariate(self.pop_params, self.indiv_params, self.observations, self.splot)
        elif self.model_type == "multivariate":
            return Multivariate(self.pop_params, self.indiv_params, self.number_param, self.observations, self.splot)
        elif self.model_type == "network":
            def compute_vtk(self, kernel_size, file_name):
                [Kxd, invKd] = make_matrices(self.dist_mat, self.ind_ctrl, kernel_size)
                signals = compute_fast_network_signal(self.pop_params, self.indiv_params, -1, Kxd, invKd)
                #write_vtk_file(self.data_mesh, signals, file_name)

            def read_matlab_data(self, path_to_matlab_data):
                mat = scipy.io.loadmat(path_to_matlab_data)
                self.dist_mat = mat["DistMat"]
                self.data_mesh = mat["data_mesh"]
                self.mesh_points = self.data_mesh["pts_left"][0][0]
                self.mesh_triangles = self.data_mesh["tri_left"][0][0]
                self.ind_ctrl = np.hstack(mat["ind_ctrl"]) - 1


    def hist(self, variable_name, nb_beans):
        plt.hist(self.pop_params[variable_name], nb_beans)
        plt.show()

    def plot_mean(self):
        return self.model.plot_mean()


    def plot_patients(self, list_of_patients_to_plot, with_obs):
        return self.model.plot_patients(list_of_patients_to_plot, with_obs)

    def hold_plot(self):
        plt.show()



    ''' UTILS '''
    #def list_pop_params(self):
    #    list_var = []
    #    for key, value in self.pop_params.iteritems():
    #        list_var.append(key)
    #    return list_var



    ''' MODEL SPECIFIC FUNCTIONS '''

