
from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')
from make_sto_basis import make_sto_basis 

name = "IKE"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_"+name+"_stochastic_0"
eta_file = "eta_used_for_read_back_"+name+"_stochastic_0"

n_sample = 50
test_node_x = [3.2375397324999998e+05, 4.0413550053600001e+05, 6.5208921350399998e+05]
test_node_y = [3.2675570484799999e+06, 3.3039835199199999e+06, 3.2880523029299998e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6374", "4509", "3879"]

dist_name = "uniform"
sto_poly_deg = 3
sto_poly_dim = 2
coefficient = [0.8, 1.2, 0.9, 1.1]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
n_modes = basis["n_modes"]

time_step = 301

mesh = Mesh(mesh_dir + mesh_file)
#mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)

eta, eta_f, u, u_f = [], [], [], []
for mode in range(n_modes):
    eta.append(Function(B))
    u.append(Function(C))
    eta_input_file = eta_file + str(mode) + ".h5"
    eta_f.append(HDF5File(mesh.mpi_comm(), input_dir + eta_input_file, "r"))
    u_input_file = u_file + str(mode) + ".h5"
    u_f.append(HDF5File(mesh.mpi_comm(), input_dir + u_input_file, "r"))

test_nodes = [[a, b] for a in test_node_x for b in test_node_y]
sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
sample_y = np.linspace(coefficient[2], coefficient[3], num = n_sample)
bin_random_u1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
bin_random_v1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
bin_random_eta1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])


for i, q0 in enumerate(sample_x):
    for j, q1 in enumerate(sample_y):
        
        print "start: q0 = " + str(round(q0, 2)) + "; q1 = " + str(round(q1, 2)) 
        orth_list = [orth[mode](q0, q1) for mode in range(n_modes)]

        for k in range(time_step):

            dataset_u = "WaterVelocity/vector_%d"%k
            dataset_eta = "SurfaceElevation/vector_%d"%k

            for mode in range(n_modes):
                u_f[mode].read(u[mode], dataset_u)
                eta_f[mode].read(eta[mode], dataset_eta)

            for l in range(len(test_nodes)):
                u1_list = [u[p](test_nodes[l][0], test_nodes[l][1])[0] for p in range(n_modes)]
                v1_list = [u[p](test_nodes[l][0], test_nodes[l][1])[1] for p in range(n_modes)]
                eta1_list = [eta[p](test_nodes[l][0], test_nodes[l][1]) for p in range(n_modes)]

                bin_random_u1[i, j, k, l] = np.dot(orth_list, u1_list)
                bin_random_v1[i, j, k, l] = np.dot(orth_list, v1_list)
                bin_random_eta1[i, j, k, l] = np.dot(orth_list, eta1_list)


np.save(output_dir + name + "_bin_random_eta1_surrogate_all_points_order_1", bin_random_eta1)
np.save(output_dir + name + "_bin_random_u1_surrogate_all_points_order_1", bin_random_u1)
np.save(output_dir + name + "_bin_random_v1_surrogate_all_points_order_1", bin_random_v1)


