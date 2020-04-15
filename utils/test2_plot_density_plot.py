

from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')
from make_sto_basis import make_sto_basis 
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
      'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick

name = "test2"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "inlet_adh_sswm_finer.xml"
u_file = "u_used_for_read_back_"+name.lower()+"_stochastic_"
eta_file = "eta_used_for_read_back_"+name.lower()+"_stochastic_"
emin = -0.24
emax = 0.24
umin = -0.2
umax = 0.2
vmin = -0.01
vmax = 0.01
form = ".png"

n_sample = 5000
bar_height = 50

test_node_x = [25.0]
test_node_y = [25.0]
test_node_str =["1/4"]
test_node_str_node_number = ["1/4"]

dist_name = "uniform"
sto_poly_deg = 3
sto_poly_dim = 2 
coefficient = [0.8, 1.2, 1.0, 2.0]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
JointCDF = basis.get("joint_cdf")
n_modes = basis["n_modes"]

time_step = 100
step_size = 0.5

#mesh = Mesh(mesh_dir + mesh_file)
mesh = RectangleMesh(Point(0, 0), Point(100, 50), 20, 10)
#mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
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
#    import pdb;pdb.set_trace()


    for i, item in enumerate(test_nodes): 

        print "    test_nodes: " + str(i) + " : x = " + str(round(item[0], 2)) + "; y = " + str(round(item[1], 2)) 

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

        #import pdb;pdb.set_trace()
        # plot
        for j, field in enumerate([u_random_output, v_random_output, eta_random_output]):
 #           import pdb;pdb.set_trace()
            try:
                plt.figure(figsize=[11.5, 7])
                plt.ylim(0, bar_height)
                plt.ylabel('Probability density function')
                if j == 0:
                    sns.distplot(field,hist_kws={"label": str(n_sample)+" samples"})
                    plt.xlim(umin,umax)
                    plt.xticks(np.linspace(umin, umax, 6))
                    plt.xlabel('x-direction velocity:m/s')
                    plt.title('Pdf of x-direction water velocity at time = ' + str(k * step_size) + ' sec')
                    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
                    plt.legend()
                    plt.savefig(output_dir + 'u_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+form)
                elif j == 1:
                    sns.distplot(field,hist_kws={"label": str(n_sample)+" samples"})
                    plt.xlim(vmin, vmax)
                    plt.xticks(np.linspace(vmin, vmax, 6))
                    plt.xlabel('y-direction velocity:m/s')
                    plt.title('Pdf of y-direction water velocity at time = ' + str(k * step_size) + ' sec')
                    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
                    plt.legend()
                    plt.savefig(output_dir + 'v_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+form)
                else:
                    sns.distplot(field,hist_kws={"label": str(n_sample)+" samples"})
                    plt.xlim(emin, emax)
                    plt.xticks(np.linspace(emin, emax, 6))
                    plt.xlabel('Surface Elevation:m')
                    plt.title('Pdf of surface elevation at time = ' + str(k * step_size) + ' sec')
                    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
                    plt.legend()
                    plt.savefig(output_dir + 'eta_display_density_plot_testnode_' + str(item[0]) + '_' + str(item[1]) +'_timestep_'+str(k)+form)
                plt.close()
            except:
                print "this field formulate a singular matrix. skip this time step at present."
                pass
