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
case = 'HARVEY'
name = 'v'
model_file = case + '_bin_random_'+name+'1_surrogate_all_points_order_1.npy'
truth_file = case + '_bin_random_'+name+'1_true_all_points_order_1.npy'
a = np.load(loc + case + dir + model_file)
b = np.load(loc + case + dir + truth_file)

x1 = 0.8 
x2 = 1.2
tick_num = 6
n_sample = 50

time_step = 1160 
dt = 447
nnode = 21 

for i in range(time_step):
    fig, ax = plt.subplots(figsize = [11.5, 7])
    X = np.linspace(x1, x2, n_sample)
    vmin = min(np.min(a[:, i, nnode]), np.min(b[:, i, nnode]))
    vmax = max(np.max(a[:, i, nnode]), np.max(b[:, i, nnode]))
    if vmin < 0:
        vmin = 1.2 * vmin
    else:
        vmin = 0.8 * vmin
    if vmax < 0:
        vmax = 0.8 * vmax
    else:
        vmax = 1.2 * vmax
    print vmin, vmax
    ax.plot(X[:], a[:, i, nnode], '-', label = 'SSWM', color = 'mediumseagreen', linewidth = 4)
    ax.plot(X[:], b[:, i, nnode], '*', label = "truth", color = 'darkorange', markersize = 12)
    ax.legend(loc = "best")
    plt.xlabel(r'$\xi_1$')
    if name == 'eta':
        plt.ylabel('Surface Elevation:m')
    elif name == 'u':
        plt.ylabel('x-direction Water Velocity:m/s')
    elif name == 'v':
        plt.ylabel('y-direction Water Velocity:m/s')
    plt.xlim([x1, x2])
    plt.ylim([vmin, vmax])
    plt.xticks(np.linspace(x1, x2, tick_num))
    plt.yticks(np.linspace(vmin, vmax, tick_num))
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.title('time:'+ '{:.3f}'.format(i * dt / 3600.0 / 24.0)+' days')
    plt.savefig(loc + case + out_dir + name + '_display_surrogate_comparison_'+str(i)+'.pdf', bbox_inches='tight')
    plt.close()



