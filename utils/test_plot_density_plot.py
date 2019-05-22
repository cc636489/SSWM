

from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')
from make_sto_basis import make_sto_basis 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick

name = "IKE"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_gulf_winds_ike_stochastic_"
eta_file = "eta_used_for_read_back_gulf_winds_ike_stochastic_"

n_sample = 1000
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
JointCDF = basis.get("joint_cdf")
n_modes = basis["n_modes"]

time_step = 500

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
