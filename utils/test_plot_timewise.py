import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mtick


loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/'
dir = '_bins/'
out_dir = '_results/'
case = 'HARVEY'
name = 'eta'
model_file = case + '_bin_random_'+name+'1_surrogate_all_points_order_1.npy'
truth_file = case + '_bin_random_'+name+'1_true_all_points_order_1.npy'
time_file = 'time_stamp_at_every_time_step.npy'
a = np.load(loc + case + dir + model_file)
b = np.load(loc + case + dir + truth_file)
t = np.load(loc + case + dir + time_file)
t = [t[i]/3600.0/24.0 for i in range(len(t))]
#import pdb;pdb.set_trace()
coefficient = [0.8, 1.2]
n_sample = 50
tick_num = 10
nnode = 21

if len(coefficient) == 4:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)
    sample_y = np.linspace(coefficient[2], coefficient[3], num = n_sample)
elif len(coefficient) == 2:
    sample_x = np.linspace(coefficient[0], coefficient[1], num = n_sample)

for i, q0 in enumerate(sample_x):
    #for j, q1 in enumerate(sample_y):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 6, forward = True)
        if len(coefficient) == 4:
            ax.plot(t[:], a[i, j, :, nnode], '-', label = 'model solution')
            ax.plot(t[:], b[i, j, :, nnode], 'x', label = 'true solution')
            vmin = min(np.min(a[i, j, :, nnode]), np.min(b[i, j, :, nnode]))
            vmax = max(np.max(a[i, j, :, nnode]), np.max(b[i, j, :, nnode]))
        elif len(coefficient) == 2:
            ax.plot(t[:], a[i, :, nnode], '-', label = 'model solution')
            ax.plot(t[:], b[i, :, nnode], 'x', label = 'true solution')
            vmin = min(np.min(a[i, :, nnode]), np.min(b[i, :, nnode]))
            vmax = max(np.max(a[i, :, nnode]), np.max(b[i, :, nnode]))
        ax.legend(loc = "upper right")
        if vmin < 0:
            vmin = 1.2 * vmin
        else:
            vmin = 0.8 * vmin
        if vmax < 0:
            vmax = 0.8 * vmax
        else:
            vmax = 1.2 * vmax
        print vmin, vmax
        plt.xlim([min(t), max(t)])
        plt.ylim([vmin, vmax])
        plt.xticks(np.linspace(min(t), max(t), tick_num))
        #import pdb;pdb.set_trace()
        plt.yticks(np.linspace(vmin, vmax, tick_num))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
        plt.xlabel('time:day')
        if name == 'eta':
            plt.ylabel('Surface Elevation:m')
        elif name == 'u':
            plt.ylabel('x-direction Water Velocity:m/s')
        elif name == 'v':
            plt.ylabel('y-direction Water Velocity:m/s')
        if len(coefficient) == 4:
            plt.title("Comparison over time at random point (" + str(round(q0, 4)) + " ," + str(round(q1, 4)) + ")")
            plt.savefig(loc + case + out_dir + name + '_display_timeseries_comparison_' + str(round(q0, 4)) + '_' + str(round(q1, 4)) + '.png')
        elif len(coefficient) == 2:
            plt.title("Comparison over time at random point q0=" + str(round(q0, 4)))
            plt.savefig(loc + case + out_dir + name + '_display_timeseries_comparison_' + str(round(q0, 4)) + '.png')

        plt.close()
 
