import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import vtk

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


###################################
### Export to paraview
###################################

def convert_signal(signal, data_mesh):
    indices = data_mesh["indices"][0][0]
    Q = data_mesh["Q"][0][0]

    output_signal = np.zeros((data_mesh["pts_left"][0][0].shape[0], 1))

    for i, idx in enumerate(indices):
        OK = np.where(Q == idx)
        output_signal[OK] = signal[i]

    return output_signal

def write_vtk_file(data_mesh, signal, file_name):

    ### Initialize parameters
    mesh_points = data_mesh["pts_left"][0][0]
    mesh_polygones = data_mesh["tri_left"][0][0]
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

    output_signal = convert_signal(signal[0], data_mesh)

    cm = plt.get_cmap('jet')
    print output_signal.shape, mesh_polygones.shape
    for i, col in enumerate(output_signal):
        colors.InsertNextTuple3(int(col[0]*255/np.max(output_signal)),0,0)


    polydata.GetPointData().SetScalars(colors)
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
