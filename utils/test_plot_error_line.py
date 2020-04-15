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
x1 = 0.8
x2 = 1.2
x3 = 0.9
x4 = 1.1
num = 10
m = (x2-x1)/num
n = (x4-x3)/num

index = [3, 7]

scalex = 0.9999
scaley = 1.0001

x = np.linspace(x1, x2, num)
y = np.linspace(x3, x4, num)

for i in range(a.shape[2]):
    for inx in index:

        plt.figure(figsize=[10, 7])
        vmin = min(min(a[inx, :, i]), min(b[inx, :, i]))
        vmax = max(max(a[inx, :, i]), max(b[inx, :, i]))
        plt.plot(x, a[inx, :, i], '-^', x, b[inx, :, i], '*')
        plt.legend(['model solution', 'true solution'])
        plt.xlabel(r'$\xi_1$')
        plt.ylabel(name)

        if min(x)>0 and max(x)>0:
            plt.xlim([min(x)*scalex, max(x)*scaley])
            plt.xticks(np.linspace(min(x)*scalex, max(x)*scaley, 11))
        elif min(x)<0 and max(x)>0:
            plt.xlim([min(x)*scaley, max(x)*scaley])
            plt.xticks(np.linspace(min(x) * scaley, max(x) *scaley, 11))
        elif min(x)<0 and max(x)<0:
            plt.xlim([min(x) * scaley, max(x) * scalex])
            plt.xticks(np.linspace(min(x) * scaley, max(x) * scalex, 11))

        if vmin>0 and vmax>0:
            plt.ylim([vmin*scalex, vmax*scaley])
            plt.yticks(np.linspace(vmin*scalex, vmax*scaley, 11))
        elif vmin<0 and vmax>0:
            plt.ylim([vmin*scaley, vmax * scaley])
            plt.yticks(np.linspace(vmin*scaley, vmax*scaley, 11))
        elif vmin<0 and vmax<0:
            plt.ylim([vmin*scaley, vmax * scalex])
            plt.yticks(np.linspace(vmin*scaley, vmax*scalex, 11))

        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))
        plt.title('time=' + str(i) + ' sec; ' + r'$\xi_2$=' + '{:0.3f}'.format(y[inx]))
        # plt.show()
        plt.savefig(loc+dir+name+'_time_'+str(i)+'_index_' + str(inx) + '_xi-2_refine_color_range'+'.png')

        plt.figure(figsize=[10, 7])
        vmin = min(min(a[:, inx, i]), min(b[:, inx, i]))
        vmax = max(max(a[:, inx, i]), max(b[:, inx, i]))
        plt.plot(y, a[:, inx, i], '-^', y, b[:, inx, i], '*')
        plt.legend(['model solution', 'true solution'])
        plt.xlabel(r'$\xi_2$')
        plt.ylabel(name)

        if min(y) > 0 and max(y) > 0:
            plt.xlim([min(y) * scalex, max(y) * scaley])
            plt.xticks(np.linspace(min(y) * scalex, max(y) * scaley, 11))
        elif min(y) < 0 and max(y) > 0:
            plt.xlim([min(y) * scaley, max(y) * scaley])
            plt.xticks(np.linspace(min(y) * scaley, max(y) * scaley, 11))
        elif min(y) < 0 and max(y) < 0:
            plt.xlim([min(y) * scaley, max(y) * scalex])
            plt.xticks(np.linspace(min(y) * scaley, max(y) * scalex, 11))

        if vmin > 0 and vmax > 0:
            plt.ylim([vmin * scalex, vmax * scaley])
            plt.yticks(np.linspace(vmin * scalex, vmax * scaley, 11))
        elif vmin < 0 and vmax > 0:
            plt.ylim([vmin * scaley, vmax * scaley])
            plt.yticks(np.linspace(vmin * scaley, vmax * scaley, 11))
        elif vmin < 0 and vmax < 0:
            plt.ylim([vmin * scaley, vmax * scalex])
            plt.yticks(np.linspace(vmin * scaley, vmax * scalex, 11))

        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4e'))
        plt.title('time=' + str(i) + ' sec; ' + r'$\xi_1$=' + '{:0.3f}'.format(x[inx]))
        # plt.show()
        plt.savefig(loc + dir + name + '_time_' + str(i) + '_index_' + str(inx) + '_xi-1_refine_color_range' + '.png')



