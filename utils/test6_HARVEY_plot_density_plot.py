

from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')
from make_sto_basis import make_sto_basis 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick

name = "HARVEY"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_gulf_winds_harvey_stochastic_4.0_2.0_"
eta_file = "eta_used_for_read_back_gulf_winds_harvey_stochastic_4.0_2.0_"

n_sample = 1000

test_node_x = [3.2402186690199998e+05, 4.6935285107099998e+05, 6.6736137265599996e+05, 3.2856664622000000e+05, 3.2529322481400002e+05, 3.1858816293200001e+05, 3.1292030319900002e+05, 3.3593948564500001e+05, 3.0007843549300003e+05, 3.1927603514400002e+05, 2.9937036285099998e+05, 2.8432231212800002e+05, 3.4697452542100003e+05, 4.1839055220600002e+05, 5.8269936381600006e+05, 6.4208428986400005e+05, 6.7038572594100004e+05, 7.7422587433200004e+05, 8.0185519664700003e+05, 8.4953476710299996e+05, 9.2267523367999995e+05, 2.5763542973900001e+05, 1.5587698318400001e+05, 2.1535995006000000e+05, 1.3629092187399999e+05, 1.1259518368700000e+05, 1.5817665928699999e+05, 1.3834004919600001e+05, 8.6808859704999995e+04, 6.3161031436899997e+04]
test_node_y = [3.2672211939200000e+06, 3.3140610491200001e+06, 3.2889479893000000e+06, 3.2655972475200002e+06, 3.2666293017500001e+06, 3.2677337144300002e+06, 3.2832635089699998e+06, 3.2879437652500002e+06, 3.2967217365700002e+06, 3.3048165327400002e+06, 3.2577057229499999e+06, 3.2374947815700001e+06, 3.2788214789999998e+06, 3.2966514931999999e+06, 3.2805300291300002e+06, 3.2774865211399999e+06, 3.2702405453200000e+06, 3.2443050481599998e+06, 3.3604420403800001e+06, 3.2448157875399999e+06, 3.2986607205500002e+06, 3.2118147654100000e+06, 3.1535162250700002e+06, 3.1849672177400002e+06, 3.1397782485500001e+06, 3.1266345019399999e+06, 3.1217132363399998e+06, 3.1076120203700000e+06, 3.0972237951000002e+06, 3.0523672295300001e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6373", "4508", "3878", "6027", "6316", "6630", "7547", "7727", "7795", "7855", "6380", "5441", "5170", "4359", "3515", "3701", "3524", "2811", "5925", "2111", "1560", "4660", "4436", "4307", "4553", "4742", "4139", "4294", "4646", "4641"]

dist_name = "uniform"
sto_poly_deg = 1
sto_poly_dim = 1
coefficient = [0.8, 1.2]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
JointCDF = basis.get("joint_cdf")
n_modes = basis["n_modes"]

time_step = 1160

mesh = Mesh(mesh_dir + mesh_file)
#mesh = RectangleMesh(Point(0, 0), Point(100, 50), 20, 10)
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

for k in range(time_step):

    print "start: timestep = " + str(k)

    dataset_u = "WaterVelocity/vector_%d"%k
    dataset_eta = "SurfaceElevation/vector_%d"%k

    for mode in range(n_modes):
        u_f[mode].read(u[mode], dataset_u)
        eta_f[mode].read(eta[mode], dataset_eta)

    for i, item in enumerate(test_nodes): 

        print "    test_nodes: " + str(i) + " : x = " + str(round(item[0], 2)) + "; y = " + str(round(item[1], 2)) 

        u1_list = [u[p](item[0], item[1])[0] for p in range(n_modes)]
        v1_list = [u[p](item[0], item[1])[1] for p in range(n_modes)]
        eta1_list = [eta[p](item[0], item[1]) for p in range(n_modes)]

        u_random_output, v_random_output, eta_random_output = [], [], []
        for m in range(n_sample):
            if len(coefficient) == 4:
                orth_list = [orth[mode](samples[0][m], samples[1][m]) for mode in range(n_modes)]
            elif len(coefficient) == 2:
                orth_list = [orth[mode](samples[m]) for mode in range(n_modes)]
            u_random_output.append(np.dot(orth_list, u1_list))
            v_random_output.append(np.dot(orth_list, v1_list))
            eta_random_output.append(np.dot(orth_list, eta1_list))

        # plot
        for j, field in enumerate([u_random_output, v_random_output, eta_random_output]):
            #import pdb;pdb.set_trace()
            try:
                plt.figure(figsize=[10, 6])
                sns.distplot(field, bins=30)
                plt.xticks(np.linspace(min(field), max(field), 10))
                plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
                plt.ylabel('Probability density function')
                if j == 0:
                    plt.xlabel('x-direction velocity')
                    plt.title('Pdf of x-direction water velocity at node ('+ str(item[0]) + ', ' + str(item[1]) +'); time step = ' + str(k) + ';')
                    plt.savefig(output_dir + 'u_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+'.png')
                elif j == 1:
                    plt.xlabel('y-direction velocity')
                    plt.title('Pdf of y-direction water velocity at node ('+ str(item[0]) + ', ' + str(item[1]) +'); time step = ' + str(k) + ';')
                    plt.savefig(output_dir + 'v_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+'.png')
                else:
                    plt.xlabel('Surface Elevation')
                    plt.title('Pdf of surface elevation at node ('+ str(item[0]) + ', ' + str(item[1]) +'); time step = ' + str(k) + ';')
                    plt.savefig(output_dir + 'eta_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+'.png')
                plt.close()
            except:
                print "this field formulate a singular matrix. skip this time step at present."
                pass
