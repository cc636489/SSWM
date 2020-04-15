
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
                      'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
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


#file1 = "8771341_Galveston_Bay_Entrance_North_Jetty_29_21.4N_94_43.5W.csv"
#file2 = "8768094_Calcasieu_Pass_29_46.1N_93_20.6W.csv"
#file3 = "8764227_LAWMA_Amerada_Pass_29_27N_91_20.3W.csv"


test_node_x = [3.2402186690199998e+05, 4.6935285107099998e+05, 6.6736137265599996e+05, 3.2856664622000000e+05, 3.2529322481400002e+05, 3.1858816293200001e+05, 3.1292030319900002e+05, 3.3593948564500001e+05, 3.0007843549300003e+05, 3.1927603514400002e+05, 2.9937036285099998e+05, 2.8432231212800002e+05, 3.4697452542100003e+05, 4.1839055220600002e+05, 5.8269936381600006e+05, 6.4208428986400005e+05, 6.7038572594100004e+05, 7.7422587433200004e+05, 8.0185519664700003e+05, 8.4953476710299996e+05, 9.2267523367999995e+05, 2.5763542973900001e+05, 1.5587698318400001e+05, 2.1535995006000000e+05, 1.3629092187399999e+05, 1.1259518368700000e+05, 1.5817665928699999e+05, 1.3834004919600001e+05, 8.6808859704999995e+04, 6.3161031436899997e+04]
test_node_y = [3.2672211939200000e+06, 3.3140610491200001e+06, 3.2889479893000000e+06, 3.2655972475200002e+06, 3.2666293017500001e+06, 3.2677337144300002e+06, 3.2832635089699998e+06, 3.2879437652500002e+06, 3.2967217365700002e+06, 3.3048165327400002e+06, 3.2577057229499999e+06, 3.2374947815700001e+06, 3.2788214789999998e+06, 3.2966514931999999e+06, 3.2805300291300002e+06, 3.2774865211399999e+06, 3.2702405453200000e+06, 3.2443050481599998e+06, 3.3604420403800001e+06, 3.2448157875399999e+06, 3.2986607205500002e+06, 3.2118147654100000e+06, 3.1535162250700002e+06, 3.1849672177400002e+06, 3.1397782485500001e+06, 3.1266345019399999e+06, 3.1217132363399998e+06, 3.1076120203700000e+06, 3.0972237951000002e+06, 3.0523672295300001e+06]
test_node_str = ["8771341", "8768094", "8764227"]
test_node_str_node_number = ["6373", "4508", "3878", "6027", "6316", "6630", "7547", "7727", "7795", "7855", "6380", "5441", "5170", "4359", "3515", "3701", "3524", "2811", "5925", "2111", "1560", "4660", "4436", "4307", "4553", "4742", "4139", "4294", "4646", "4641"]


time_step = 1160
time = np.load(time_dir + time_file)
time[:] = time[:]/3600.0/24.0

# read in ADCIRC data.
tmp = []
for item in test_node_str_node_number:
    tmp.append(pd.read_csv(time_dir + "ADCIRC_" + item + "0.csv"))

# read in observation time and predicted.
time_obs = tmp[0]["Time"].tolist()
time_obs = [time_obs[i]/3600.0/24.0 for i in range(len(time_obs))]
tmp_obs = []
for i in range(len(test_node_str_node_number)):
    tmp_obs.append(tmp[i]["elev"].tolist())

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
        plt.figure(figsize=[11.5, 7])

        cc1 = field_mean[:, i] - field_var[:, i]
        cc2 = field_mean[:, i] + field_var[:, i]
        if k == 0: 
            plt.plot(time, field_mean[:, i], '-', label='SSWM:'+r'$\mu\pm\sigma$', alpha=0.2, color='#1B2ACC')
            plt.plot(time_obs[::10], tmp_obs[i][::10], "*", label='ADCIRC', color="darkorange")
            ymin = min(min(cc1), min(tmp_obs[i]))
            ymax = max(max(cc2), max(tmp_obs[i]))
        else:
            plt.plot(time, field_mean[:, i], '-', alpha=0.2, color='#1B2ACC')
            ymin = min(cc1)
            ymax = max(cc2)
        plt.legend()
        plt.fill_between(time, field_mean[:, i] - field_var[:, i],  field_mean[:, i] + field_var[:, i], alpha=0.2,
                         edgecolor='mediumseagreen', facecolor='mediumseagreen', linewidth=0.2, linestyle='dashdot', antialiased=True)
        plt.xlim([min(time), max(time)])
        plt.xticks(np.linspace(min(time), max(time), 6))
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
            plt.ylabel('Surface Elevation: m')
            plt.title('Comparison of surface elevation at node ' + test_node_str_node_number[i])
            plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            plt.savefig(output_dir + 'eta_validation_ADCIRC_line_mean_variance_' + test_node_str_node_number[i][:7] +'.pdf')
        elif k == 1:
            plt.ylabel('x-direction Water Velocity: m/s')
            plt.title('Comparison of x-direction water velocity at node ' + test_node_str_node_number[i])
            plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            plt.savefig(output_dir + 'u_validation_ADCIRC_line_mean_variance_' + test_node_str_node_number[i][:7] +'.pdf')
        else:
            plt.ylabel('y-direction Water Velocity: m/s')
            plt.title('Comparison of y-direction water velocity at node ' + test_node_str_node_number[i])
            plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
            plt.savefig(output_dir + 'v_validation_ADCIRC_line_mean_variance_' + test_node_str_node_number[i][:7] +'.pdf')
        #plt.legend(['mean', r'mean$\pm\sigma$'])
        plt.close()
        print "Done: " + str(i) + "test_node."



