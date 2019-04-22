
import numpy as np
import matplotlib.pyplot as plt
from fenics import *


name = "test2"
#input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/500_steps_060/"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+\
             "_results/"
time_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+\
           "_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
#mesh_file = "inlet_adh_sswm_finer.xml"
#mean_u_file = "u_used_for_read_back_" + name + "_stochastic_finer_mesh_7_1_060_00.h5"
#mean_eta_file = "eta_used_for_read_back_" + name + "_stochastic_finer_mesh_7_1_060_00.h5"
mean_u_file = "u_used_for_read_back_" + name + "_stochastic_00.h5"
mean_eta_file = "eta_used_for_read_back_" + name + "_stochastic_00.h5"
time_file = "time_stamp_at_every_time_step.npy"
variance_eta_file = "eta_used_for_read_back_variance_in_domain.h5"
variance_u_file = "u_used_for_read_back_variance_in_domain.h5"
variance_v_file = "v_used_for_read_back_variance_in_domain.h5"

test_node_x = [25.0, 50.0, 75.0]
test_node_y = [25.0]

time_step = 100
time = np.load(time_dir + time_file)
#time[:] = time[:]/3600/24
time = time[:-1]

#mesh = Mesh(mesh_dir + mesh_file)
mesh = RectangleMesh(Point(0, 0), Point(100, 50), 20, 10)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
D = FunctionSpace(mesh, "CG", 2)
mean_u = Function(C)
mean_eta = Function(B)
var_u = Function(D)
var_v = Function(D)
var_eta = Function(B)
var_u_new = Function(B)
var_v_new = Function(B)

test_nodes = [[a, b] for a in test_node_x for b in test_node_y]
plot_eta_mean = np.zeros([time_step, len(test_nodes)])
plot_u_mean = np.zeros([time_step, len(test_nodes)])
plot_v_mean = np.zeros([time_step, len(test_nodes)])
plot_eta_var = np.zeros([time_step, len(test_nodes)])
plot_u_var = np.zeros([time_step, len(test_nodes)])
plot_v_var = np.zeros([time_step, len(test_nodes)])

mean_u_f = HDF5File(mesh.mpi_comm(), input_dir + mean_u_file, "r")
mean_eta_f = HDF5File(mesh.mpi_comm(), input_dir + mean_eta_file, "r")
variance_eta_f = HDF5File(mesh.mpi_comm(), output_dir + variance_eta_file, "r")
variance_u_f = HDF5File(mesh.mpi_comm(), output_dir + variance_u_file, "r")
variance_v_f = HDF5File(mesh.mpi_comm(), output_dir + variance_v_file, "r")

for k in range(time_step):

    dataset_mean_u = "WaterVelocity/vector_%d"%k
    dataset_mean_eta = "SurfaceElevation/vector_%d"%k
    dataset_variance_u = "UVelocity/vector_%d"%k
    dataset_variance_v = "VVelocity/vector_%d"%k
    dataset_variance_eta = "SurfaceElevation/vector_%d"%k

    mean_u_f.read(mean_u, dataset_mean_u)
    mean_eta_f.read(mean_eta, dataset_mean_eta)
    variance_eta_f.read(var_eta, dataset_variance_eta)
    variance_u_f.read(var_u, dataset_variance_u)
    variance_v_f.read(var_v, dataset_variance_v)

    # if np.isnan(var_eta.vector().array()).any():
    #     print "sha bi"
    #     import pdb;pdb.set_trace()
    # if np.isnan(var_u.vector().array()).any():
    #     print "sha bi"
    #     import pdb;pdb.set_trace()
    # if np.isnan(var_v.vector().array()).any():
    #     print "sha bi"
    #     import pdb;pdb.set_trace()

    var_u_new.interpolate(var_u)
    var_v_new.interpolate(var_v)
    var_u_new.vector()[abs(var_u_new.vector().array()) < 1e-10] = 0
    var_v_new.vector()[abs(var_v_new.vector().array()) < 1e-10] = 0
    var_eta.vector()[abs(var_eta.vector().array()) < 1e-10] = 0

    # if not all([var_eta.vector().array()[i] >= 0 for i in range(len(var_eta.vector().array()))]):
    #     print "sha bi"
    #     import pdb;pdb.set_trace()
    # if not all([var_u.vector().array()[i] >= 0 for i in range(len(var_u.vector().array()))]):
    #     print "sha bi"
    #     import pdb;pdb.set_trace()
    # if not all([var_v.vector().array()[i] >= 0 for i in range(len(var_v.vector().array()))]):
    #     print "sha bi"
    #     import pdb;pdb.set_trace()

    # if not all([var_u_new.vector().array()[i] >= 0 for i in range(len(var_u_new.vector().array()))]):
    #     print "sha bi"
    #     import pdb;pdb.set_trace()
    # if not all([var_v_new.vector().array()[i] >= 0 for i in range(len(var_v_new.vector().array()))]):
    #     print "sha bi"
    #     import pdb;pdb.set_trace()

    # def lt0(x):
    #     if x < 0:
    #         return True
    #     return False
    # vresult = filter(lt0, var_v_new.vector().array())
    # uresult = filter(lt0, var_u_new.vector().array())
    # import pdb;pdb.set_trace()


    for i, item in enumerate(test_nodes):
        plot_eta_mean[k, i] = mean_eta(item[0], item[1])
        plot_u_mean[k, i] = mean_u(item[0], item[1])[0]
        plot_v_mean[k, i] = mean_u(item[0], item[1])[1]
        plot_eta_var[k, i] = abs(var_eta(item[0], item[1])) ** 0.5
        plot_u_var[k, i] = abs(var_u_new(item[0], item[1])) ** 0.5
        plot_v_var[k, i] = abs(var_v_new(item[0], item[1])) ** 0.5

for k, [field_mean, field_var] in enumerate([[plot_eta_mean, plot_eta_var], [plot_u_mean, plot_u_var], [plot_v_mean,
                                                                                                         plot_v_var]]):
    for i in range(len(test_nodes)):
        plt.figure(figsize=[7, 6])
        plt.plot(time, field_mean[:, i], '-', alpha=0.2, color='#1B2ACC')
        plt.fill_between(time, field_mean[:, i] - field_var[:, i],  field_mean[:, i] + field_var[:, i], alpha=0.2,
                         edgecolor='#1B2ACC', facecolor='#089FFF', linewidth=0.2, linestyle='dashdot', antialiased=True)
        # plt.errorbar(time, mean, yerr=field_std, fmt='-^', color='#FF9848')
        plt.xlim([min(time), max(time)])
        plt.xticks(np.linspace(min(time), max(time), 10))
        if min(field_mean[:, i] - field_var[:, i]) <= 0 and max(field_mean[:, i] + field_var[:, i]) <= 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 1.2, max(field_mean[:, i] + field_var[:, i]) * 0.8])
        elif min(field_mean[:, i] - field_var[:, i]) <= 0 and max(field_mean[:, i] + field_var[:, i]) > 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 1.2, max(field_mean[:, i] + field_var[:, i]) * 1.2])
        elif min(field_mean[:, i] - field_var[:, i]) > 0 and max(field_mean[:, i] + field_var[:, i]) > 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 0.8, max(field_mean[:, i] + field_var[:, i]) * 1.2])
        else:
            print min(field_mean[:, i] - field_var[:, i]), max(field_mean[:, i] + field_var[:, i])
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 0.8, max(field_mean[:, i] + field_var[:, i]) * 0.8])
        plt.xlabel('time:sec')
        if k == 0:
            plt.ylabel('Surface Elevation:m')
            plt.title('Mean and variance of surface elevation at physical point (' + str(test_nodes[i][0]) + ', ' +
                      str(test_nodes[i][1]) + ")" )
            plt.savefig(output_dir + 'eta_display_line_mean_variance_' + str(test_nodes[i][0]) + '_' +
                        str(test_nodes[i][1]) + '.png')
        elif k == 1:
            plt.ylabel('x-direction Water Velocity:m/s')
            plt.title('Mean and variance of x-direction water velocity at physical point (' + str(test_nodes[i][0]) +
                      ', ' + str(test_nodes[i][1]) + ")" )
            plt.savefig(output_dir + 'u_display_line_mean_variance_' + str(test_nodes[i][0]) + '_' +
                        str(test_nodes[i][1]) + '.png')
        else:
            plt.ylabel('y-direction Water Velocity:m/s')
            plt.title('Mean and variance of y-direction water velocity at physical point (' + str(test_nodes[i][0]) +
                      ', ' + str(test_nodes[i][1]) + ")" )
            plt.savefig(output_dir + 'v_display_line_mean_variance_' + str(test_nodes[i][0]) + '_' +
                        str(test_nodes[i][1]) + '.png')
        plt.legend(['mean', r'mean$\pm\sigma$'])
        plt.close()
        # plt.show()
        print "Done: " + str(i) + "test_node."



