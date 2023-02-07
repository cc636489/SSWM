
from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/Users/chenchen/gloria/7-SSWM-github/src/')
from make_sto_basis import make_sto_basis 
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
      'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick

name = "test3"
#input_dir = "/Users/chenchen/gloria/4-NS-results-and-tests/test2_res/convergence_inprobability_test_result_dimension_2/"
input_dir = "/Users/chenchen/gloria/test_fenics/stochastic/test3_stochastic_hump/"
output_dir = "/Users/chenchen/gloria/test_fenics/stochastic/test3_stochastic_hump/"
mesh_dir = "/Users/chenchen/gloria/7-SSWM-github/input/"
st_order = 3
u_file = "u_used_for_read_back_"+name.lower()+"_stochastic_"
eta_file = "eta_used_for_read_back_"+name.lower()+"_stochastic_"

n_sample = 5000

test_node_x = [500.0]
test_node_y = [100.0]
test_node_str =["1/4"]
test_node_str_node_number = ["1/4"]

dist_name = "uniform"
sto_poly_deg = st_order
sto_poly_dim = 2
coefficient = [0.8, 1.2, 0.9, 1.1]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
JointCDF = basis.get("joint_cdf")
n_modes = basis["n_modes"]

time_step = 155 + 1
step_size = 0.5

#mesh = Mesh(mesh_dir + mesh_file)
#mesh = RectangleMesh(Point(0, 0), Point(100, 50), 20, 10)
mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)

eta, eta_f, u, u_f = [], [], [], []
for mode in range(n_modes):
    eta.append(Function(B))
    u.append(Function(C))
    eta_input_file = eta_file + '{:02d}'.format(mode) + ".h5"
    eta_f.append(HDF5File(mesh.mpi_comm(), input_dir + eta_input_file, "r"))
    u_input_file = u_file + '{:02d}'.format(mode) + ".h5"
    u_f.append(HDF5File(mesh.mpi_comm(), input_dir + u_input_file, "r"))

#test_nodes = [[a, b] for a in test_node_x for b in test_node_y]
test_nodes = [[test_node_x[i], test_node_y[i]] for i in range(len(test_node_x))]
samples = JointCDF.sample(n_sample) 

for k in range(time_step-1, time_step):

    print("start: timestep = " + str(k))

    dataset_u = "WaterVelocity/vector_%d"%k
    dataset_eta = "SurfaceElevation/vector_%d"%k

    for mode in range(n_modes):
        u_f[mode].read(u[mode], dataset_u)
        eta_f[mode].read(eta[mode], dataset_eta)
#    import pdb;pdb.set_trace()


    for i, item in enumerate(test_nodes): 

        print("    test_nodes: " + str(i) + " : x = " + str(round(item[0], 2)) + "; y = " + str(round(item[1], 2)))

        u1_list = [u[p](item[0], item[1])[0] for p in range(n_modes)]
        v1_list = [u[p](item[0], item[1])[1] for p in range(n_modes)]
        eta1_list = [eta[p](item[0], item[1]) for p in range(n_modes)]

        #import pdb;pdb.set_trace()
        u_random_output, v_random_output, eta_random_output = [], [], []
        for m in range(n_sample):
            if len(coefficient) == 4:
                orth_list = [orth[mode](samples[0][m], samples[1][m]) for mode in range(n_modes)]
            elif len(coefficient) == 2:
                orth_list = [orth[mode](samples[m]) for mode in range(n_modes)]
            u_random_output.append(np.dot(orth_list, u1_list))
            v_random_output.append(np.dot(orth_list, v1_list))
            eta_random_output.append(np.dot(orth_list, eta1_list))


        # save to file
        np.save(output_dir + "binrandom_surrogate_u_samples_order_"+str(item[0])+"_"+str(item[1])+"_"+str(st_order)+"_NSample_"+str(n_sample), u_random_output)
        np.save(output_dir + "binrandom_surrogate_v_samples_order_"+str(item[0])+"_"+str(item[1])+"_"+str(st_order)+"_NSample_"+str(n_sample), v_random_output)
        np.save(output_dir + "binrandom_surrogate_eta_samples_order_"+str(item[0])+"_"+str(item[1])+"_"+str(st_order)+"_NSample_"+str(n_sample), eta_random_output)


       
