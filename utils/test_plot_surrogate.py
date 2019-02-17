import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mtick


loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
dir = 'test4_res/test4_res_perturb_sto2_fx_1e-6q0q1/'
name = 'v'
model_file = 'binrandom_'+name+'1_order_2.npy'
truth_file = 'binrandom_'+name+'1_true_order_2.npy'
a = np.load(loc+dir+model_file)
b = np.load(loc+dir+truth_file)
# vmin = min(np.min(a), np.min(b))
# vmax = max(np.max(a), np.max(b))
x1 = 0.8
x2 = 1.2
x3 = 0.9
x4 = 1.1
num = 10.
m = (x2-x1)/num
n = (x4-x3)/num

for i in range(a.shape[2]):
    fig, ax = plt.subplots()
    X = np.linspace(x1, x2, num)
    Y = np.linspace(x3, x4, num)
    X, Y = np.meshgrid(X, Y)
    vmin = np.min(a[:, :, i])
    vmax = np.max(a[:, :, i])
    # plt.pcolor(X, Y, a[:, :, i], cmap=cm.jet, vmin=vmin, vmax=vmax)
    plt.pcolor(X, Y, a[:, :, i], cmap=cm.jet, vmin=vmin, vmax=vmax)

    X = np.reshape(X, [1, -1])
    Y = np.reshape(Y, [1, -1])
    A = np.reshape(b[:, :, i], [1, -1])
    plt.scatter(X, Y, c=A, cmap=cm.jet, marker='o', s=30, linewidth=0.5, edgecolors='black', vmin=vmin, vmax=vmax)
    plt.grid()
    plt.xlabel(r'$\xi_1$')
    plt.ylabel(r'$\xi_2$')
    plt.xlim([x1, x2])
    plt.ylim([x3, x4])
    plt.xticks(np.linspace(x1, x2, num))
    plt.yticks(np.linspace(x3, x4, num))
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    plt.title('time:'+str(i)+' sec')
    cbar = plt.colorbar(format='%.5e')
    if name == 'eta':
        cbar.ax.set_title(r'$\eta_{model}, \eta_{truth}$')
    elif name == 'u':
        cbar.ax.set_title(r'$u_{model}, u_{truth}$')
    elif name == 'v':
        cbar.ax.set_title(r'$v_{model}, v_{truth}$')
    # plt.show()
    plt.savefig(loc+dir+name+'_display_relax_color_range'+str(i)+'.png')



