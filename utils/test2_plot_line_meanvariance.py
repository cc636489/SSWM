
import numpy as np
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

import matplotlib.pyplot as plt
from fenics import *
from matplotlib.ticker import FormatStrFormatter


name = "test2"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results_largest_variance/"
output_dir = input_dir 

#name = "HARVEY"
#input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
#output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
time_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
#mesh_file = "inlet_adh_sswm_finer.xml"
#mean_u_file = "u_used_for_read_back_" + name + "_stochastic_0.2_0.75_0.13les_00.h5"
#mean_eta_file = "eta_used_for_read_back_" + name + "_stochastic_0.2_0.75_0.13les_00.h5"
mesh_file = "Gulf_wind.xml"
#mean_u_file = "u_used_for_read_back_" + name.lower() + "_stochastic_00.h5"
#mean_eta_file = "eta_used_for_read_back_" + name.lower() + "_stochastic_00.h5"
mean_u_file = "u_used_for_read_back_" + name.lower() + "_stochastic_largest_variance_00.h5"
mean_eta_file = "eta_used_for_read_back_" + name.lower() + "_stochastic_largest_variance_00.h5"
#mean_u_file = "u_used_for_read_back_gulf_winds_" + name.lower() + "_stochastic_4.0_2.0_00.h5"
#mean_eta_file = "eta_used_for_read_back_gulf_winds_" + name.lower() + "_stochastic_4.0_2.0_00.h5"
time_file = "time_stamp_at_every_time_step.npy"
variance_eta_file = "eta_used_for_read_back_variance_in_domain.h5"
variance_u_file = "u_used_for_read_back_variance_in_domain.h5"
variance_v_file = "v_used_for_read_back_variance_in_domain.h5"

#test_node_x = [1635473.5225233955, 1653973.1434681825, 1684252.1148920991, 1708037.3418211113, 1705117.6737536343, 1706149.6253981739, 1742746.1544508543]
#test_node_y = [3342589.6136079477, 3278394.675424758, 3313831.7655779063, 3281734.2964863107, 3268004.7432332593, 3262809.777137509, 3070224.9625879396]
#test_node_str = ["8761927, New Canal Station", "8764227, LAWMA Amerada Pass", "8768094, Calcasieu Pass", "8771013, Galveston Bay-Eagle Point", "8771341, Galveston Bay Entrance-North Jetty", "8771450, Galveston Pier", "8775870, Bob Hall Pier"]

#test_node_x = [3.2402186690199998e+05, 4.6935285107099998e+05, 6.6736137265599996e+05, 3.2856664622000000e+05, 3.2529322481400002e+05, 3.1858816293200001e+05, 3.1292030319900002e+05, 3.3593948564500001e+05, 3.0007843549300003e+05, 3.1927603514400002e+05, 2.9937036285099998e+05, 2.8432231212800002e+05, 3.4697452542100003e+05, 4.1839055220600002e+05, 5.8269936381600006e+05, 6.4208428986400005e+05, 6.7038572594100004e+05, 7.7422587433200004e+05, 8.0185519664700003e+05, 8.4953476710299996e+05, 9.2267523367999995e+05, 2.5763542973900001e+05, 1.5587698318400001e+05, 2.1535995006000000e+05, 1.3629092187399999e+05, 1.1259518368700000e+05, 1.5817665928699999e+05, 1.3834004919600001e+05, 8.6808859704999995e+04, 6.3161031436899997e+04]
#test_node_y = [3.2672211939200000e+06, 3.3140610491200001e+06, 3.2889479893000000e+06, 3.2655972475200002e+06, 3.2666293017500001e+06, 3.2677337144300002e+06, 3.2832635089699998e+06, 3.2879437652500002e+06, 3.2967217365700002e+06, 3.3048165327400002e+06, 3.2577057229499999e+06, 3.2374947815700001e+06, 3.2788214789999998e+06, 3.2966514931999999e+06, 3.2805300291300002e+06, 3.2774865211399999e+06, 3.2702405453200000e+06, 3.2443050481599998e+06, 3.3604420403800001e+06, 3.2448157875399999e+06, 3.2986607205500002e+06, 3.2118147654100000e+06, 3.1535162250700002e+06, 3.1849672177400002e+06, 3.1397782485500001e+06, 3.1266345019399999e+06, 3.1217132363399998e+06, 3.1076120203700000e+06, 3.0972237951000002e+06, 3.0523672295300001e+06]
#test_node_str = ["8771341", "8768094", "8764227"]
#test_node_str_node_number = ["6373", "4508", "3878", "6027", "6316", "6630", "7547", "7727", "7795", "7855", "6380", "5441", "5170", "4359", "3515", "3701", "3524", "2811", "5925", "2111", "1560", "4660", "4436", "4307", "4553", "4742", "4139", "4294", "4646", "4641"]

test_node_x = [25.0, 50.0, 75.0]
test_node_y = [25.0, 25.0, 25.0]
test_node_str = ["first quadrant", "middle", "third quadrant"]
test_node_str_node_number = ["1","2","3"]

time_step = 100
time = np.load(time_dir + time_file)
#time[:] = time[:]/3600/24
#time = time[:-1]

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
#import pdb; pdb.set_trace()
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
        fig, ax = plt.subplots(figsize=[11.5, 8])
        plt.plot(time, field_mean[:, i], '-', alpha=0.2, color='#1B2ACC')
        plt.fill_between(time, field_mean[:, i] - field_var[:, i],  field_mean[:, i] + field_var[:, i], alpha=0.2,
                         edgecolor='#1B2ACC', facecolor='#089FFF', linewidth=0.2, linestyle='dashdot', antialiased=True)
        # plt.errorbar(time, mean, yerr=field_std, fmt='-^', color='#FF9848')
        #plt.xlim([min(time), max(time)])
        plt.xlim([0, 50])
        #plt.xticks(np.linspace(min(time), max(time), 10))
        plt.xticks(np.linspace(0, 50, 5))
        if min(field_mean[:, i] - field_var[:, i]) <= 0 and max(field_mean[:, i] + field_var[:, i]) <= 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 1.2, max(field_mean[:, i] + field_var[:, i]) * 0.8])
        elif min(field_mean[:, i] - field_var[:, i]) <= 0 and max(field_mean[:, i] + field_var[:, i]) > 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 1.2, max(field_mean[:, i] + field_var[:, i]) * 1.2])
        elif min(field_mean[:, i] - field_var[:, i]) > 0 and max(field_mean[:, i] + field_var[:, i]) > 0:
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 0.8, max(field_mean[:, i] + field_var[:, i]) * 1.2])
        else:
            print min(field_mean[:, i] - field_var[:, i]), max(field_mean[:, i] + field_var[:, i])
            plt.ylim([min(field_mean[:, i] - field_var[:, i]) * 0.8, max(field_mean[:, i] + field_var[:, i]) * 0.8])
       # plt.xlabel('time:day')
        plt.xlabel('time:sec')
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        if k == 0:
            plt.ylabel('Surface Elevation:m')
            plt.title('Deviation of elevation at point (' + str(test_node_x[i]) + ', ' +
                      str(test_node_y[i]) + ")" )
            plt.savefig(output_dir + 'eta_display_line_mean_largest_variance_' + str(test_node_x[i]) + '_' +
                        str(test_node_y[i]) + '.pdf')
            #plt.title('Mean and variance of surface elevation at node ' + test_node_str_node_number[i])
            #plt.savefig(output_dir + 'eta_display_line_mean_variance_at_node_' + test_node_str_node_number[i] +'.png')
        elif k == 1:
            plt.ylabel('x-direction Water Velocity:m/s')
            plt.title('Deviation of u velocity at point (' + str(test_node_x[i]) +
                      ', ' + str(test_node_y[i]) + ")" )
            plt.savefig(output_dir + 'u_display_line_mean_largest_variance_' + str(test_node_x[i]) + '_' +
                        str(test_node_y[i]) + '.pdf')
            #plt.title('Mean and variance of x-direction water velocity at station ' + test_node_str_node_number[i])
            #plt.savefig(output_dir + 'u_display_line_mean_variance_at_node_' + test_node_str_node_number[i] +'.png')
        else:
            plt.ylabel('y-direction Water Velocity:m/s')
            plt.title('Deviation of v velocity at point (' + str(test_node_x[i]) +
                      ', ' + str(test_node_y[i]) + ")" )
            plt.savefig(output_dir + 'v_display_line_mean_largest_variance_' + str(test_node_x[i]) + '_' +
                        str(test_node_y[i]) + '.pdf')
            #plt.title('Mean and variance of y-direction water velocity at station ' + test_node_str_node_number[i])
            #plt.savefig(output_dir + 'v_display_line_mean_variance_at_node_' + test_node_str_node_number[i] +'.png')
        plt.legend(['mean', r'mean$\pm\sigma$'])
        plt.close()
        # plt.show()
        print "Done: " + str(i) + "test_node."



