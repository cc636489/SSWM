import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mtick


loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/'
dir = '_bins/'
out_dir = '_results/'
case = 'test4'
name = 'v'
model_file = case + '_bin_random_'+name+'1_surrogate_all_points_order_1.npy'
truth_file = case + '_bin_random_'+name+'1_true_all_points_order_1.npy'
time_file = 'time_stamp_at_every_time_step.npy'
a = np.load(loc + case + dir + model_file)
b = np.load(loc + case + dir + truth_file)
t = np.load(loc + case + dir + time_file)
t = [t[i]/3600.0/24.0 for i in range(len(t))]
#import pdb;pdb.set_trace()
x1 = 0.9
x2 = 1.1
x3 = 0.9
x4 = 1.1
n_sample = 10
tick_num = 10
sample_x = np.linspace(x1, x2, num = n_sample)
sample_y = np.linspace(x3, x4, num = n_sample)

for i, q0 in enumerate(sample_x):
    for j, q1 in enumerate(sample_y):
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 6, forward = True)
        ax.plot(t[:], a[i, j, :, 0], '-', label = 'model solution')
        ax.plot(t[:], b[i, j, :, 0], 'x', label = 'true solution')
        ax.legend(loc = "upper right")
        vmin = min(np.min(a[i, j, :, 0]), np.min(b[i, j, :, 0]))
        vmax = max(np.max(a[i, j, :, 0]), np.max(b[i, j, :, 0]))
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
        plt.title("Comparison over time at random point (" + str(round(q0, 2)) + " ," + str(round(q1, 2)) + ")")

        #plt.show()
        plt.close()
        plt.savefig(loc + case + out_dir + name + '_display_timeseries_comparison_' + str(round(q0, 2)) + '_' + str(round(q1, 2)) + '.png')
 
