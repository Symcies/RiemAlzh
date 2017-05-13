### Imports

import numpy as np
import scipy.io
import vtk

###################################
### Reading function
###################################
def read_population_parameters(file_name):
    params = {}
    f = open(file_name, 'r')

    while True:
        # Read line, continue if empty and break if last one
        line = f.readline()
        if not line:
            break
        if line == "":
            continue

        # Preprocessing
        line = line.lower()

        # Extract parameters
        list_elements = line.split(' ')
        param_name = list_elements[0]
        param_values = []
        for i in range(1, len(list_elements)):
            try:
                param_values.append(float(list_elements[i]))
            except:
                pass

        # Add to existing dict
        params[param_name] = param_values

    return params

def add_individual_parameters(line, number_of_params):
    individual_parameters = {}

    list_elements = line.split(' ')
    idx_start = 0

    for tpl in number_of_params:
        idx_end = idx_start + tpl[1]
        individual_parameters[tpl[0].rstrip()] = [float(list_elements[i].rstrip()) for i in range(idx_start, idx_end)]
        idx_start = idx_end

    return individual_parameters

def read_individual_parameters(file_name):
    params = []
    labels = []
    number_per_label = []
    individuals = []

    f = open(file_name, 'r')

    line_count = -1

    while True:
        line_count += 1
        # Read line, continue if empty and break if last one
        line = f.readline()
        if not line:
            break
        if line == "":
            continue

        # Read line, continue if empty, break if last one
        if (line_count == 0):
            labels = line.split(' ')
            continue

        if (line_count == 1):
            number_per_label = line.split(' ')

            for idx, name in enumerate(labels):
                params.append((name, int(number_per_label[idx])))

            continue

        # Add individuals
        individuals.append(add_individual_parameters(line, params))

    return individuals


###################################
### Reading matlab files
###################################

mat = scipy.io.loadmat('CorticalSurface.mat')
dist_mat = mat["DistMat"]
data_mesh = mat["data_mesh"]
mesh_points = data_mesh["pts_left"][0][0]
mesh_triangles = data_mesh["tri_left"][0][0]
ind_ctrl = np.hstack(mat["ind_ctrl"])
ind_ctrl = ind_ctrl - 1 # There is a difference within the indexation between matlab and python

def make_matrices(dist_mat, ind_ctrl, bandwidth):
    N_v = dist_mat.shape[0]
    N_c = len(ind_ctrl)

    K = np.exp( - dist_mat**2 / bandwidth**2)
    Kd = K[np.ix_(ind_ctrl, ind_ctrl)]
    invKd = np.linalg.inv(Kd)

    Kxd = np.zeros((N_v, N_c))

    for i in range(0, N_v):
        for j in range(0, N_c):
            Kxd[i, j] = K[i, ind_ctrl[j]]

    return [Kxd, invKd]

[Kxd, invKd] = make_matrices(dist_mat, ind_ctrl, 10)

pop_file = "Work_pop.txt"
indiv_file = "Work_indiv.txt"

pop_params = read_population_parameters(pop_file)
indiv_params = read_individual_parameters(indiv_file)

###################################
### Computing the signals for the networks
###################################

def compute_fast_network_signal(pop_params, indiv_params, indiv_number, Kxd, invKd):

    #delta = compute_interpolation(Kxd, invKd, pop_params["delta"])
    #nu = compute_interpolation(Kxd, invKd, pop_params["nu"])
    delta = pop_params["delta"]
    nu = pop_params["nu"]

    thickness = pop_params["thickness"][0]


    time_points = np.linspace(60,100,20)
    space_shift = [0 for i in range(0, len(delta))]
    if(indiv_number > -1):
        time_points = np.exp(indiv_params["ksi"][0]) * (time_points - indiv_params["tau"][0])
        space_shift = indiv_params["w"]

    signals = []
    for t in time_points:
        signal = np.divide(space_shift, thickness*np.exp(delta)) + delta + np.divide(nu, thickness * t)
        signal = thickness * np.exp(signal)
        signals.append(signal)

    return signals


signals = compute_fast_network_signal(pop_params, indiv_params, -1, Kxd, invKd)


###################################
### Export to paraview
###################################

def convert_signal(signal, data_mesh):
    indices = data_mesh["indices"]
    Q = data_mesh["Q"]

    output_signal = arrayofsizeNumberOfPoints

    for i, idx in enumerate(indices):
        ### indices where Q == idx
        ### output_signal[indices_where_Q == idx] = signal(i)

    return output_signal

def write_vtk_file(mesh_points, mesh_polygones, signal, file_name):

    ### Initialize parameters
    nb_points, dimension = mesh_points.shape
    nb_polygones, polygone_dimensions = mesh_polygones.shape
    # TODO : Assert : dimension should be 3
    # TODO : Assert : POints not empty
    # TODO : Assert : Polygone dimension
    # TODO : Assert : polygone not empty

    ### Initialize the output file
    points = vtk.vtkPoints()
    vectices = vtk.vtkCellArray()

    ### Write the points coordinate
    for pt in mesh_points:
        points.InsertNextPoint((pt[0], pt[1], pt[2]))

    ### Write the triangles
    for pl in mesh_polygones:
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
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")

    for i in range(0, len(mesh_polygones)):
        colors.InsertNextTuple3(255,0,0)

    polydata.GetCellData().SetScalars(colors)
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


write_vtk_file(mesh_points, mesh_triangles, 1, "mesh.vtp")
