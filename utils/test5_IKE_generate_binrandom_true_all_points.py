
from fenics import *
import numpy as np
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')

name = "IKE"
inout_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_truth/"
out_dir = inout_dir[:-6] + "bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_gulf_winds_ike_deterministic_"
eta_file = "eta_used_for_read_back_gulf_winds_ike_deterministic_"

n_sample = 50
test_node_x = [3.2375397324999998e+05, 4.0413550053600001e+05, 6.5208921350399998e+05, 3.2900153713800001e+05, 3.2557148308400001e+05, 3.1881131950699998e+05, 3.1460159233700001e+05, 3.3496542305799999e+05, 3.0139283829799999e+05, 3.1282448336200003e+05, 2.9570583575099998e+05, 2.8491055709199997e+05, 3.4851441155299998e+05, 4.2815074394299998e+05, 5.7302018531800003e+05, 6.3071504403200001e+05, 6.6204347742899996e+05, 7.8584841874100000e+05, 7.9657833684899996e+05, 8.6587532801900001e+05, 9.3000227841499995e+05, 2.6392368590400001e+05, 1.6298658246000001e+05]
test_node_y = [3.2675570484799999e+06, 3.3039835199199999e+06, 3.2880523029299998e+06, 3.2657850455399998e+06, 3.2670047864700002e+06, 3.2679612539499998e+06, 3.2836322031399999e+06, 3.2890447269899999e+06, 3.2982392603799999e+06, 3.3043017858099998e+06, 3.2589011959699998e+06, 3.2381142812700002e+06, 3.2793550391199999e+06, 3.2975682191800000e+06, 3.2830367486999999e+06, 3.2823205112999999e+06, 3.2709672468599998e+06, 3.2448157875399999e+06, 3.3663629658800000e+06, 3.2474602108100001e+06, 3.3006520252700001e+06, 3.2122789727400001e+06, 3.1579897588000000e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6374", "4509", "3879", "6028", "6317", "6631", "7548", "7728", "7796", "7856", "6381", "5442", "5171", "4360", "3516", "3702", "3525", "2812", "5926", "2112", "1561", "4661", "4437"]

coefficient = [0.8, 1.2]

time_step = 500

mesh = Mesh(mesh_dir + mesh_file)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
eta = Function(B)
u = Function(C)

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
    #for j, q1 in enumerate(sample_y):
        
        #print "start: q0 = " + str(round(q0, 2)) + "; q1 = " + str(round(q1, 2)) 
        print "start: q0 = " + str(round(q0, 4))

        #u_input_file = u_file + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_00.h5"
        u_input_file = u_file + str(round(q0, 4)) + "_00.h5"
        u_f = HDF5File(mesh.mpi_comm(), inout_dir + u_input_file, "r")
        for k in range(time_step):
            dataset = "WaterVelocity/vector_%d"%k
            attr = u_f.attributes(dataset)
            u_f.read(u, dataset)
            for l in range(len(test_nodes)):
                if len(coefficient) == 4:
                    bin_random_u1[i, j, k, l] = u(test_nodes[l][0], test_nodes[l][1])[0]
                    bin_random_v1[i, j, k, l] = u(test_nodes[l][0], test_nodes[l][1])[1]
                elif len(coefficient) == 2:
                    bin_random_u1[i, k, l] = u(test_nodes[l][0], test_nodes[l][1])[0] 
                    bin_random_v1[i, k, l] = u(test_nodes[l][0], test_nodes[l][1])[1] 

        #eta_input_file = eta_file + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_00.h5"
        eta_input_file = eta_file + str(round(q0, 4)) + "_00.h5"
        eta_f = HDF5File(mesh.mpi_comm(), inout_dir + eta_input_file, "r")
        for k in range(time_step):
            dataset = "SurfaceElevation/vector_%d"%k
            attr = eta_f.attributes(dataset)
            eta_f.read(eta, dataset)
            for l in range(len(test_nodes)):
                if len(coefficient) == 4:
                    bin_random_eta1[i, j, k, l] = eta(test_nodes[l][0], test_nodes[l][1])
                elif len(coefficient) == 2:
                    bin_random_eta1[i, k, l] = eta(test_nodes[l][0], test_nodes[l][1])
        
np.save(out_dir + name + "_bin_random_eta1_true_all_points_order_1", bin_random_eta1)
np.save(out_dir + name + "_bin_random_u1_true_all_points_order_1", bin_random_u1)
np.save(out_dir + name + "_bin_random_v1_true_all_points_order_1", bin_random_v1)


