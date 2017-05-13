from visu import visu

pop_file = "Work_pop.txt"
indiv_file = "Work_indiv.txt"


network = visu("fast-network", "data/Work_pop.txt", "data/Work_indiv.txt", 1)
network.hist("delta", 50)
#network.read_matlab_data("data/CorticalSurface.mat")
#network.compute_vtk(16, "mesh2.vtp")
