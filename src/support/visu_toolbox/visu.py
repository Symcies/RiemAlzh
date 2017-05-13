import numpy as np
import matplotlib.pyplot as plt
import scipy.io



from reader import read_individual_parameters, read_population_parameters
from fast_network import make_matrices, compute_fast_network_signal, write_vtk_file

class visu:
    '''
    Documentation
    '''



    def __init__(self, type, population_params, individual_params, observations):
        self.type = type
        self.pop_params = read_population_parameters(population_params)
        self.indiv_params = read_individual_parameters(individual_params)
        #self.observations = read_observations(observations)

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
        return 0

    def list_pop_params(self):
        list_var = []
        for key, value in self.pop_params.iteritems():
            list_var.append(key)
        return list_var


    def plot_patients(self, list_of_patients_to_plot):
        return 0


    def compute_vtk(self, kernel_size, file_name):
        [Kxd, invKd] = make_matrices(self.dist_mat, self.ind_ctrl, kernel_size)
        signals = compute_fast_network_signal(self.pop_params, self.indiv_params, -1, Kxd, invKd)
        write_vtk_file(self.data_mesh, signals, file_name)
