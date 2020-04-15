
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fenics import *


name = "HARVEY"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
time_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "Gulf_wind.xml"
mean_u_file = "u_used_for_read_back_gulf_winds_harvey_stochastic_4.0_2.0_00.h5"
mean_eta_file = "eta_used_for_read_back_gulf_winds_harvey_stochastic_4.0_2.0_00.h5"
time_file = "time_stamp_at_every_time_step.npy"
variance_eta_file = "eta_used_for_read_back_variance_in_domain.h5"
variance_u_file = "u_used_for_read_back_variance_in_domain.h5"
variance_v_file = "v_used_for_read_back_variance_in_domain.h5"


test_node_x = [3.2375397324999998e+05, 4.0413550053600001e+05, 6.5208921350399998e+05, 1.5732570298100001e+05, 7.8602831890000001e+04, 6.2312613404500000e+04] 
test_node_y = [3.2675570484799999e+06, 3.3039835199199999e+06, 3.2880523029299998e+06, 3.1645656953099999e+06, 3.0968750273400000e+06, 3.0732804932200001e+06]
test_node_str = ["8771341", "8768094", "8764227", "8773767", "8775870", "8775241"]
test_node_str_node_number = ["6373", "4508", "3878", "4557", "4739", "4820"]
file1 = "CO-OPS_8771341_wl.csv"
file2 = "CO-OPS_8768094_wl.csv"
file3 = "CO-OPS_8764227_wl.csv"
file4 = "CO-OPS_8773767_wl.csv"
file5 = "CO-OPS_8775870_wl.csv"
file6 = "CO-OPS_8775241_wl.csv"

time_step = 1160
time = np.load(time_dir + time_file)
time[:] = time[:]/3600.0/24.0

# read in observation time.
tmp1 = pd.read_csv(time_dir + file1)
tmp2 = pd.read_csv(time_dir + file2)
tmp3 = pd.read_csv(time_dir + file3)
tmp4 = pd.read_csv(time_dir + file4)
tmp5 = pd.read_csv(time_dir + file5)
tmp6 = pd.read_csv(time_dir + file6)

# read in observation time and predicted.
time_obs = tmp1["time"].tolist()
time_obs = [time_obs[i]/3600.0/24.0 for i in range(len(time_obs))]
tmp1_obs = tmp1["verified"].tolist()
tmp2_obs = tmp2["verified"].tolist()
tmp3_obs = tmp3["verified"].tolist()
tmp4_obs = tmp4["verified"].tolist()
tmp5_obs = tmp5["verified"].tolist()
tmp6_obs = tmp6["verified"].tolist()
#import pdb; pdb.set_trace()

tmp1_obs = [0.0 if x == '-' else float(x) for x in tmp1_obs]
tmp2_obs = [0.0 if x == '-' else float(x) for x in tmp2_obs]
tmp3_obs = [0.0 if x == '-' else float(x) for x in tmp3_obs]
tmp4_obs = [0.0 if x == '-' else float(x) for x in tmp4_obs]
tmp5_obs = [0.0 if x == '-' else float(x) for x in tmp5_obs]
tmp6_obs = [0.0 if x == '-' else float(x) for x in tmp6_obs]
#import pdb;pdb.set_trace()

mesh = Mesh(mesh_dir + mesh_file)
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

# check node inside domain or not.
remove_node = []
for i in xrange(len(test_node_x)):
    try:
        mean_eta(test_node_x[i], test_node_y[i])
    except:
        print "this node i=" + str(i) + ", " + test_node_str[i] + "is outside the domain, delete it."
        remove_node.append(i)

if len(remove_node) == len(test_node_x):
    test_nodes = [[test_node_x[i], test_node_y[i]] for i in range(len(test_node_x))]
    mean_u.set_allow_extrapolation(True)
    mean_eta.set_allow_extrapolation(True)
    var_u.set_allow_extrapolation(True)
    var_v.set_allow_extrapolation(True)
    var_eta.set_allow_extrapolation(True)
    var_u_new.set_allow_extrapolation(True)
    var_v_new.set_allow_extrapolation(True)
else:
    test_nodes = [[test_node_x[i], test_node_y[i]] for i in range(len(test_node_x)) if i not in remove_node]
    test_node_str = [test_node_str[i] for i in range(len(test_node_str)) if i not in remove_node]

#test_nodes = [[a, b] for a in test_node_x for b in test_node_y]
import pdb; pdb.set_trace()
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

    var_u_new.interpolate(var_u)
    var_v_new.interpolate(var_v)
    var_u_new.vector()[abs(var_u_new.vector().array()) < 1e-10] = 0
    var_v_new.vector()[abs(var_v_new.vector().array()) < 1e-10] = 0
    var_eta.vector()[abs(var_eta.vector().array()) < 1e-10] = 0

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

        cc1 = field_mean[:, i] - field_var[:, i]
        cc2 = field_mean[:, i] + field_var[:, i]
        if k == 0 and i == 0:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp1_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp1_obs))
            ymax = max(max(cc2), max(tmp1_obs))
        elif k == 0 and i == 1:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp2_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp2_obs))
            ymax = max(max(cc2), max(tmp2_obs))
        elif k == 0 and i == 2:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp3_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp3_obs))
            ymax = max(max(cc2), max(tmp3_obs))
        elif k == 0 and i == 3:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp4_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp4_obs))
            ymax = max(max(cc2), max(tmp4_obs))
        elif k == 0 and i == 4:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp5_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp5_obs))
            ymax = max(max(cc2), max(tmp5_obs))
        elif k == 0 and i == 5:
            plt.plot(time, field_mean[:, i], '-', label='surrogate model', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp6_obs[::10], "*", label='observation', color="orange")
            ymin = min(min(cc1), min(tmp6_obs))
            ymax = max(max(cc2), max(tmp6_obs))
        else:
            plt.plot(time, field_mean[:, i], '-', alpha=0.2, color='#1B2ACC')
            ymin = min(cc1)
            ymax = max(cc2)

        plt.legend(loc="upper left")
        plt.fill_between(time, field_mean[:, i] - field_var[:, i],  field_mean[:, i] + field_var[:, i], alpha=0.2,
                         edgecolor='#1B2ACC', facecolor='#089FFF', linewidth=0.2, linestyle='dashdot', antialiased=True)
        plt.xlim([min(time), max(time)])
        plt.xticks(np.linspace(min(time), max(time), 10))
        if ymin <= 0 and ymax <= 0:
            plt.ylim([ymin * 1.2, ymax * 0.8])
        elif ymin <= 0 and ymax > 0:
            plt.ylim([ymin * 1.2, ymax * 1.2])
        elif ymin > 0 and ymax > 0:
            plt.ylim([ymin * 0.8, ymax * 1.2])
        else:
            plt.ylim([ymin * 0.8, ymax * 0.8])
        print ymin, ymax
        plt.xlabel('time:day')
        if k == 0:
            plt.ylabel('Surface Elevation:m')
            plt.title('Mean and variance of surface elevation at station ' + test_node_str[i])
            plt.savefig(output_dir + 'eta_validation_line_mean_variance_' + test_node_str[i][:14] +'.png')
        elif k == 1:
            plt.ylabel('x-direction Water Velocity:m/s')
            plt.title('Mean and variance of x-direction water velocity at station ' + test_node_str[i])
            plt.savefig(output_dir + 'u_validation_line_mean_variance_' + test_node_str[i][:14] +'.png')
        else:
            plt.ylabel('y-direction Water Velocity:m/s')
            plt.title('Mean and variance of y-direction water velocity at station ' + test_node_str[i])
            plt.savefig(output_dir + 'v_validation_line_mean_variance_' + test_node_str[i][:14] +'.png')
        plt.legend(['mean', r'mean$\pm\sigma$'])
        plt.close()
        # plt.show()
        print "Done: " + str(i) + "test_node."



