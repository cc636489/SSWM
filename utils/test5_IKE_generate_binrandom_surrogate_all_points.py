
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
u_file = "u_used_for_read_back_gulf_winds_ike_stochastic_"
eta_file = "eta_used_for_read_back_gulf_winds_ike_stochastic_"

n_sample = 50 
test_node_x = [3.2375397324999998e+05, 4.0413550053600001e+05, 6.5208921350399998e+05, 3.2900153713800001e+05, 3.2557148308400001e+05, 3.1881131950699998e+05, 3.1460159233700001e+05, 3.3496542305799999e+05, 3.0139283829799999e+05, 3.1282448336200003e+05, 2.9570583575099998e+05, 2.8491055709199997e+05, 3.4851441155299998e+05, 4.2815074394299998e+05, 5.7302018531800003e+05, 6.3071504403200001e+05, 6.6204347742899996e+05, 7.8584841874100000e+05, 7.9657833684899996e+05, 8.6587532801900001e+05, 9.3000227841499995e+05, 2.6392368590400001e+05, 1.6298658246000001e+05]
test_node_y = [3.2675570484799999e+06, 3.3039835199199999e+06, 3.2880523029299998e+06, 3.2657850455399998e+06, 3.2670047864700002e+06, 3.2679612539499998e+06, 3.2836322031399999e+06, 3.2890447269899999e+06, 3.2982392603799999e+06, 3.3043017858099998e+06, 3.2589011959699998e+06, 3.2381142812700002e+06, 3.2793550391199999e+06, 3.2975682191800000e+06, 3.2830367486999999e+06, 3.2823205112999999e+06, 3.2709672468599998e+06, 3.2448157875399999e+06, 3.3663629658800000e+06, 3.2474602108100001e+06, 3.3006520252700001e+06, 3.2122789727400001e+06, 3.1579897588000000e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6374", "4509", "3879", "6028", "6317", "6631", "7548", "7728", "7796", "7856", "6381", "5442", "5171", "4360", "3516", "3702", "3525", "2812", "5926", "2112", "1561", "4661", "4437"]


dist_name = "uniform"
sto_poly_deg = 1
sto_poly_dim = 1
coefficient = [0.8, 1.2]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
n_modes = basis["n_modes"]

time_step = 500

mesh = Mesh(mesh_dir + mesh_file)
# mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
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


test_nodes = [[test_node_x[i], test_node_y[i]] for i in range(len(test_node_x))]
if len(coefficient) == 4:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
    sample_y = np.linspace(coefficient[2], coefficient[3], num = n_sample)
    bin_random_u1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
    bin_random_v1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
    bin_random_eta1 = np.zeros([n_sample, n_sample, time_step, len(test_nodes)])
elif len(coefficient) == 2:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
    bin_random_u1 = np.zeros([n_sample, time_step, len(test_nodes)])
    bin_random_v1 = np.zeros([n_sample, time_step, len(test_nodes)])
    bin_random_eta1 = np.zeros([n_sample, time_step, len(test_nodes)])


for i, q0 in enumerate(sample_x):
#    for j, q1 in enumerate(sample_y):
        
        #print "start: q0 = " + str(round(q0, 2)) + "; q1 = " + str(round(q1, 2)) 
        print "start: q0 = " + str(round(q0, 2))

        #orth_list = [orth[mode](q0, q1) for mode in range(n_modes)]
        orth_list = [orth[mode](q0) for mode in range(n_modes)]

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

                if len(coefficient) == 4:
                    bin_random_u1[i, j, k, l] = np.dot(orth_list, u1_list)
                    bin_random_v1[i, j, k, l] = np.dot(orth_list, v1_list)
                    bin_random_eta1[i, j, k, l] = np.dot(orth_list, eta1_list)
                elif len(coefficient) == 2:
                    bin_random_u1[i, k, l] = np.dot(orth_list, u1_list)
                    bin_random_v1[i, k, l] = np.dot(orth_list, v1_list)
                    bin_random_eta1[i, k, l] = np.dot(orth_list, eta1_list)


np.save(output_dir + name + "_bin_random_eta1_surrogate_all_points_order_1", bin_random_eta1)
np.save(output_dir + name + "_bin_random_u1_surrogate_all_points_order_1", bin_random_u1)
np.save(output_dir + name + "_bin_random_v1_surrogate_all_points_order_1", bin_random_v1)


