import numpy as np
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal', 'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mtick


loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/'
dir = '_bins/'
out_dir = '_results/'
case = 'test4'
name = 'u'
model_file = case + '_bin_random_'+name+'1_surrogate_all_points_order_1.npy'
truth_file = case + '_bin_random_'+name+'1_true_all_points_order_1.npy'
time_file = 'time_stamp_at_every_time_step.npy'
a = np.load(loc + case + dir + model_file)
b = np.load(loc + case + dir + truth_file)
t = np.load(loc + case + dir + time_file)
t = [t[i]/3600.0/24.0 for i in range(len(t))]
#t = [t[i] for i in range(len(t))]
#import pdb;pdb.set_trace()
#coefficient = [0.8, 1.2, 0.9, 1.1]
coefficient = [1.0, 2.0]
n_sample = 50
tick_num = 6
nnode = 2 

if len(coefficient) == 4:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
    sample_y = np.linspace(coefficient[2], coefficient[3], num = n_sample)
elif len(coefficient) == 2:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)

for i, q0 in enumerate(sample_x):
    #for j, q1 in enumerate(sample_y):
        fig, ax = plt.subplots()
        fig.set_size_inches(11.5, 7, forward = True)
        if len(coefficient) == 4:
            ax.plot(t[:], a[i, j, :, nnode], '-', label = 'SSWM', color = 'mediumseagreen', linewidth = 4)
            ax.plot(t[:-1:2], b[i, j, ::2, nnode], '*', label = 'truth', color = 'darkorange', markersize = 12)
            vmin = min(np.min(a[i, j, :, nnode]), np.min(b[i, j, :, nnode]))
            vmax = max(np.max(a[i, j, :, nnode]), np.max(b[i, j, :, nnode]))
        elif len(coefficient) == 2:
            ax.plot(t[:], a[i, :, nnode], '-', label = 'SSWM', color = 'mediumseagreen', linewidth = 4)
            ax.plot(t[::4], b[i, ::4, nnode], '*', label = 'truth', color = 'darkorange', markersize = 12)
            vmin = min(np.min(a[i, :, nnode]), np.min(b[i, :, nnode]))
            vmax = max(np.max(a[i, :, nnode]), np.max(b[i, :, nnode]))
        ax.legend(loc = "best")
        if vmin < 0:
            vmin = 2.0 * vmin
        else:
            vmin = 0.5 * vmin
        if vmax < 0:
            vmax = 0.5 * vmax
        else:
            vmax = 2.0 * vmax
        print vmin, vmax
        plt.xlim([min(t), max(t)])
        plt.ylim([vmin, vmax])
        plt.xticks(np.linspace(min(t), max(t), tick_num))
        plt.yticks(np.linspace(vmin, vmax, tick_num))
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        #plt.xlabel('time(seconds)')
        plt.xlabel('time(days)')
        if name == 'eta':
            plt.ylabel('Surface Elevation:m')
        elif name == 'u':
            plt.ylabel('x-direction Water Velocity:m/s')
        elif name == 'v':
            plt.ylabel('y-direction Water Velocity:m/s')
        if len(coefficient) == 4:
            plt.title("Comparison at random point (" + str(round(q0, 4)) + " ," + str(round(q1, 4)) + ")")
            plt.savefig(loc + case + out_dir + name + '_display_timeseries_comparison_' + str(nnode) + '_' + str(round(q0, 4)) + '_' + str(round(q1, 4)) + '.pdf', bbox_inches = 'tight')
        elif len(coefficient) == 2:
            plt.title("Comparison at random point " + r'$\xi_1$'  + "=" + str(round(q0, 4)))
            plt.savefig(loc + case + out_dir + name + '_display_timeseries_comparison_' + str(nnode) + '_' + str(round(q0, 4)) + '.pdf', bbox_inches = 'tight')

        plt.close()
 
