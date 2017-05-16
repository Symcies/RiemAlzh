import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy
import vtk


matplotlib.use('tkagg')

from models.networks.fast_network import FastNetwork
from models.networks.meshwork     import Meshwork
from models.networks.network      import  Network

from utils.input_management import read_individual_parameters, read_population_parameters, extract_observations

class VisuNetwork:
    '''
    Documentation
    '''

    def __init__(self, type, population_params, individual_params, observations):
        self.model_type = type.lower()
        self.pop_params = read_population_parameters(population_params)
        self.number_param, self.indiv_params = read_individual_parameters(individual_params)
        self.observations = extract_observations(observations)
        self.model = self.init_model()


    def init_model(self):
        if self.model_type == "fast-network":
            return FastNetwork()
        elif self.model_type == "network":
            return Network()
        elif self.model_type == "meshwork":
            return Meshwork()

    ''' VISUALIZATION SPECIFIC FUNCTIONS '''
    def hist(self, variable_name, nb_beans):
        plt.hist(self.pop_params[variable_name], nb_beans)
        plt.show()


    def print_population_params(self):
        return [key for key, value in self.pop_params.iteritems()]


    ''' LOADING SPECIFIC FUNCTIONS '''
    def load_matlab_data(self, path_to_matlab_data):
        mat = scipy.io.loadmat(path_to_matlab_data)
        self.dist_mat = mat["DistMat"]
        self.data_mesh = mat["data_mesh"]
        self.mesh_points = self.data_mesh["pts_left"][0][0]
        self.mesh_triangles = self.data_mesh["tri_left"][0][0]
        self.ind_ctrl = np.hstack(mat["ind_ctrl"]) - 1

    def compute_interpolation_matrices(self, kernel_bandwidth):
        N_v = self.dist_mat.shape[0]
        N_c = len(self.ind_ctrl)

        K = np.exp( - self.dist_mat**2 / kernel_bandwidth**2)
        Kd = K[np.ix_(self.ind_ctrl, self.ind_ctrl)]
        invKd = np.linalg.inv(Kd)

        Kxd = np.zeros((N_v, N_c))

        for i in range(0, N_v):
            for j in range(0, N_c):
                Kxd[i, j] = K[i, self.ind_ctrl[j]]

        self.Kxd = Kxd
        self.invKd = invKd

    def compute_signal(self, time_points, indiv_number):
        signals = []
        for t in time_points:
            signals.append(self.model.compute_signal(self.pop_params, self.indiv_params, indiv_number, t))

        return signals

    def write_vtk(self, signals, file_name):
        for idx, s in enumerate(signals):
            self.write_vtk_file(s, file_name + str(idx) + '.vtp')

    ''' PRIVATE FUNCTIONS '''
    def convert_signal(self, signal, data_mesh):
        indices = data_mesh["indices"][0][0]
        Q = data_mesh["Q"][0][0]

        output_signal = np.zeros((data_mesh["pts_left"][0][0].shape[0], 1))

        for i, idx in enumerate(indices):
            OK = np.where(Q == idx)
            output_signal[OK] = signal[i]

        return output_signal

    def write_vtk_file(self, signal, file_name):

        ### Initialize parameters
        nb_points, dimension = self.mesh_points.shape
        nb_polygones, polygone_dimensions = self.mesh_triangles.shape
        # TODO : Assert : dimension should be 3
        # TODO : Assert : POints not empty
        # TODO : Assert : Polygone dimension
        # TODO : Assert : polygone not empty


        ### Initialize the output file
        points = vtk.vtkPoints()
        vectices = vtk.vtkCellArray()

        ### Write the points coordinate
        for pt in self.mesh_points:
            points.InsertNextPoint((pt[0], pt[1], pt[2]))

        ### Write the triangles
        for pl in self.mesh_triangles:
            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, pl[0]-1)
            triangle.GetPointIds().SetId(1, pl[1]-1)
            triangle.GetPointIds().SetId(2, pl[2]-1)
            vectices.InsertNextCell(triangle)

        ### Set up the mesh
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(vectices)

        ### Set up the colors
        colors = vtk.vtkDoubleArray()
        colors.SetNumberOfComponents(1)
        colors.SetName("Scalars")

        output_signal = self.convert_signal(signal, self.data_mesh)

        cm = plt.get_cmap('jet')

        for i, col in enumerate(output_signal):
            colors.InsertNextTuple(col)


        polydata.GetPointData().AddArray(colors)
        polydata.GetPointData().SetActiveVectors(colors.GetName())
        polydata.Modified()

        #### Write the file
        if vtk.VTK_MAJOR_VERSION <= 5:
            polydata.Update()

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file_name)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polydata)
        else:
            writer.SetInputData(polydata)

        writer.Write()
