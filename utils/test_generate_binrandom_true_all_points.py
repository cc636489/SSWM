
from fenics import *
import numpy as np
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')

name = "test4"
inout_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_truth/"
out_dir = inout_dir[:-6] + "bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "inlet_adh_sswm_finer.xml"
u_file = "u_used_for_read_back_"+name+"_deterministic_"
eta_file = "eta_used_for_read_back_"+name+"_deterministic_"

n_sample = 10
test_node_x = [0.0]
test_node_y = [0.0]

coefficient = [0.9, 1.1, 0.9, 1.1]

time_step = 500

mesh = Mesh(mesh_dir + mesh_file)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
eta = Function(B)
u = Function(C)

test_nodes = [[a, b] for a in test_node_x for b in test_node_y]
sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
sample_y = np.linspace(coefficient[2], coefficient[3], num = n_sample)
bin_random_u1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
bin_random_v1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
bin_random_eta1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])


for i, q0 in enumerate(sample_x):
    for j, q1 in enumerate(sample_y):
        
        print "start: q0 = " + str(round(q0, 2)) + "; q1 = " + str(round(q1, 2)) 

        u_input_file = u_file + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_00.h5"
        u_f = HDF5File(mesh.mpi_comm(), inout_dir + u_input_file, "r")
        for k in range(time_step):
            dataset = "WaterVelocity/vector_%d"%k
            attr = u_f.attributes(dataset)
            u_f.read(u, dataset)
            for l in range(len(test_nodes)):
                bin_random_u1[i, j, k, l] = u(test_nodes[l][0], test_nodes[l][1])[0]
                bin_random_v1[i, j, k, l] = u(test_nodes[l][0], test_nodes[l][1])[1]

        eta_input_file = eta_file + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_00.h5"
        eta_f = HDF5File(mesh.mpi_comm(), inout_dir + eta_input_file, "r")
        for k in range(time_step):
            dataset = "SurfaceElevation/vector_%d"%k
            attr = eta_f.attributes(dataset)
            eta_f.read(eta, dataset)
            for l in range(len(test_nodes)):
                bin_random_eta1[i, j, k, l] = eta(test_nodes[l][0], test_nodes[l][1])
        
np.save(out_dir + name + "_bin_random_eta1_true_all_points_order_1", bin_random_eta1)
np.save(out_dir + name + "_bin_random_u1_true_all_points_order_1", bin_random_u1)
np.save(out_dir + name + "_bin_random_v1_true_all_points_order_1", bin_random_v1)


