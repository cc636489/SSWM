
from fenics import *
import numpy as np
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')

name = "HARVEY"
inout_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_truth/"
out_dir = inout_dir[:-6] + "bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_gulf_winds_harvey_deterministic_"
eta_file = "eta_used_for_read_back_gulf_winds_harvey_deterministic_"

n_sample = 50
test_node_x = [3.2402186690199998e+05, 4.6935285107099998e+05, 6.6736137265599996e+05, 3.2856664622000000e+05, 3.2529322481400002e+05, 3.1858816293200001e+05, 3.1292030319900002e+05, 3.3593948564500001e+05, 3.0007843549300003e+05, 3.1927603514400002e+05, 2.9937036285099998e+05, 2.8432231212800002e+05, 3.4697452542100003e+05, 4.1839055220600002e+05, 5.8269936381600006e+05, 6.4208428986400005e+05, 6.7038572594100004e+05, 7.7422587433200004e+05, 8.0185519664700003e+05, 8.4953476710299996e+05, 9.2267523367999995e+05, 2.5763542973900001e+05, 1.5587698318400001e+05, 2.1535995006000000e+05, 1.3629092187399999e+05, 1.1259518368700000e+05, 1.5817665928699999e+05, 1.3834004919600001e+05, 8.6808859704999995e+04, 6.3161031436899997e+04]
test_node_y = [3.2672211939200000e+06, 3.3140610491200001e+06, 3.2889479893000000e+06, 3.2655972475200002e+06, 3.2666293017500001e+06, 3.2677337144300002e+06, 3.2832635089699998e+06, 3.2879437652500002e+06, 3.2967217365700002e+06, 3.3048165327400002e+06, 3.2577057229499999e+06, 3.2374947815700001e+06, 3.2788214789999998e+06, 3.2966514931999999e+06, 3.2805300291300002e+06, 3.2774865211399999e+06, 3.2702405453200000e+06, 3.2443050481599998e+06, 3.3604420403800001e+06, 3.2448157875399999e+06, 3.2986607205500002e+06, 3.2118147654100000e+06, 3.1535162250700002e+06, 3.1849672177400002e+06, 3.1397782485500001e+06, 3.1266345019399999e+06, 3.1217132363399998e+06, 3.1076120203700000e+06, 3.0972237951000002e+06, 3.0523672295300001e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6373", "4508", "3878", "6027", "6316", "6630", "7547", "7727", "7795", "7855", "6380", "5441", "5170", "4359", "3515", "3701", "3524", "2811", "5925", "2111", "1560", "4660", "4436", "4307", "4553", "4742", "4139", "4294", "4646", "4641"]


coefficient = [0.8, 1.2]

time_step = 1160

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


