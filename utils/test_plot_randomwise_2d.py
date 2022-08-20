import numpy as np
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal', 'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mtick


loc = '/Users/chenchen/gloria/4-NS-results-and-tests/regression_test_stochastic/'
dir = '_bins/'
out_dir = '_results_ErrorMap/'
case = 'test3'
name = 'u'
colorbar_format = '%.2e'
model_file = case + '_bin_random_'+name+'1_surrogate_all_points_order_1.npy'
truth_file = case + '_bin_random_'+name+'1_true_all_points_order_1.npy'
a = np.load(loc + case + dir + model_file)
b = np.load(loc + case + dir + truth_file)

x1 = 0.8
x2 = 1.2
x3 = 0.9
x4 = 1.1
num = 6
n_sample = 20

time_step = 300
delta_t = 1

for i in [155]:  # range(time_step):
    fig, ax = plt.subplots()
    X = np.linspace(x1, x2, n_sample)
    Y = np.linspace(x3, x4, n_sample)
    X, Y = np.meshgrid(X, Y)
    # vmin = min(np.min(a[:, :, i, 0]), np.min(b[:, :, i, 0]))
    # vmax = max(np.max(a[:, :, i, 0]), np.max(b[:, :, i, 0]))
    vmin = np.min(a[:, :, i, 0]-b[:, :, i, 0])
    vmax = np.max(a[:, :, i, 0]-b[:, :, i, 0])
    print(vmin)
    print(vmax)
    # plt.pcolor(X, Y, a[:, :, i, 0], cmap=cm.jet, vmin=vmin, vmax=vmax, shading='auto')
    plt.pcolor(X, Y, a[:, :, i, 0]-b[:, :, i, 0], cmap=cm.jet, vmin=vmin, vmax=vmax, shading='auto')

    # X = np.reshape(X, [1, -1])
    # Y = np.reshape(Y, [1, -1])
    # A = np.reshape(b[:, :, i, 0], [1, -1])
    # plt.scatter(X, Y, c=A, cmap=cm.jet, marker='o', s=30, linewidth=0.5, edgecolors='black', vmin=vmin, vmax=vmax)
    plt.grid()
    plt.xlabel(r'$\xi_1$')
    plt.ylabel(r'$\xi_2$')
    plt.xlim([x1, x2])
    plt.ylim([x3, x4])
    plt.xticks(np.linspace(x1, x2, num))
    plt.yticks(np.linspace(x3, x4, num))
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    plt.title('time:'+str(round(i*delta_t, 1))+' sec')
    cbar = plt.colorbar(format=colorbar_format, ticks=np.linspace(vmin, vmax, num))
    if name == 'eta':
        cbar.ax.set_ylabel(r'$\eta_{sswm}-\eta_{benchmark}$')
    elif name == 'u':
        cbar.ax.set_ylabel(r'$u_{sswm}-u_{benchmark}$')
    elif name == 'v':
        cbar.ax.set_ylabel(r'$v_{sswm}, v_{benchmark}$')
    #plt.show()
    # plt.savefig(loc + case + out_dir + name + '_display_surrogate_comparison_'+str(i)+'.pdf',bbox_inches='tight')
    plt.savefig(loc + case + out_dir + name + '_display_surrogate_errorMap_'+str(i)+'.pdf',bbox_inches='tight')
    plt.close()



