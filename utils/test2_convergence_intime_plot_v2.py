import numpy as np
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal', 'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test2_res/'
modestatus = 'convergence_intime_test_result/'
eta_l2_prefix = 'error_eta_l2_'
uv_l2_prefix = 'error_u_l2_'
eta_h1_prefix = 'error_eta_h1_'
uv_h1_prefix = 'error_u_h1_'
mesh_size = "timestepsize: 0.5s, 0.25s, 0.125s, 0.0625s, 0.03125s"

time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap = 0, 1, 1
sur_fix = 'initial_time_too_fine'

time_step = 1
base_dt = 0.05
which_dt = 10 
hsize = [1.0, 2.0, 4.0, 8.0, 16.0]

npoint = len(hsize) 
feta_l2 = inputdir + testcase + modestatus + eta_l2_prefix
fuv_l2 = inputdir + testcase + modestatus + uv_l2_prefix
feta_h1 = inputdir + testcase + modestatus + eta_h1_prefix
fuv_h1 = inputdir + testcase + modestatus + uv_h1_prefix
eta_l2, u_l2 = np.zeros([npoint, time_step]), np.zeros([npoint, time_step])
eta_h1, u_h1 = np.zeros([npoint, time_step]), np.zeros([npoint, time_step])

for l, k in enumerate(hsize):

    with open(feta_l2 + str(k), 'r') as f:
        for i, line in enumerate(f):
            if i == int(k*which_dt):
                eta_l2[l][0] = float(line)

    with open(fuv_l2 + str(k), 'r') as f:
        for i, line in enumerate(f):
            if i == int(k*which_dt):
                u_l2[l][0] = float(line)

    with open(feta_h1 + str(k), 'r') as f:
        for i, line in enumerate(f):
            if i == int(k*which_dt):
                eta_h1[l][0] = float(line)

    with open(fuv_h1 + str(k), 'r') as f:
        for i, line in enumerate(f):
            if i == int(k*which_dt):
                u_h1[l][0] = float(line)


temp_x = [np.log(base_dt / hsize[i]) for i in range(npoint)]


def average_slope(field, h):
    ll = len(field)
    su = 0.0
    for ii in range(ll-1):
        su += (np.log(field[ii]) - np.log(field[ii+1])) / (np.log(h[ii]) - np.log(h[ii + 1]))
    return su / (ll-1)


u_h1_slope = []
eta_h1_slope = []
u_l2_slope = []
eta_l2_slope = []
for ii in range(npoint-1):
    u_h1_slope.append((np.log(u_h1[ii, 0]) - np.log(u_h1[ii + 1, 0])) / (np.log(1.0/hsize[ii]) - np.log(1.0/hsize[ii + 1])))
    u_l2_slope.append((np.log(u_l2[ii, 0]) - np.log(u_l2[ii + 1, 0])) / (np.log(1.0/hsize[ii]) - np.log(1.0/hsize[ii + 1])))
    eta_h1_slope.append((np.log(eta_h1[ii, 0]) - np.log(eta_h1[ii + 1, 0])) / (np.log(1.0/hsize[ii]) - np.log(1.0/hsize[ii + 1])))
    eta_l2_slope.append((np.log(eta_l2[ii, 0]) - np.log(eta_l2[ii + 1, 0])) / (np.log(1.0/hsize[ii]) - np.log(1.0/hsize[ii + 1])))


print u_h1_slope
print u_l2_slope
print eta_h1_slope
print eta_l2_slope


# plots
plt.subplots(figsize=(11.5, 7))
for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap):
    temp_y = [np.log(u_l2[i, i_step]) for i in range(npoint)]
    plt.plot(temp_x, temp_y, "*-", markersize = 12, linewidth = 4, color = 'mediumseagreen')
plt.legend(['slope=' + str(round(average_slope(u_l2[:, i_step], hsize), 4)) for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap)], loc='upper left')
plt.grid(which='both')
plt.xlim([-4, 1])
plt.xticks(np.linspace(-4, 1, 6))
plt.ylim([-15, -7])
plt.yticks(np.linspace(-15, -7, 6))
plt.xlabel('log(h)')
plt.ylabel('log('+ r'$\Vert e_u  \Vert$'+ ')-L2norm')
plt.title('L2 convergence test for velocity')
# plt.show()
plt.savefig(fuv_l2 + 'convergence_' + sur_fix + '.pdf', bbox_inches='tight')

plt.subplots(figsize=(11.5, 7))
for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap):
    temp_y = [np.log(u_h1[i, i_step]) for i in range(npoint)]
    plt.plot(temp_x, temp_y, "*-", markersize = 12, linewidth = 4, color = 'mediumseagreen')
plt.legend(['slope=' + str(round(average_slope(u_h1[:, i_step], hsize), 4)) for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap)], loc='upper left')
plt.grid(which='both')
plt.xlim([-4, 1])
plt.xticks(np.linspace(-4, 1, 6))
plt.ylim([-11, -6])
plt.yticks(np.linspace(-11, -6, 6))
plt.xlabel('log(h)')
plt.ylabel('log('+ r'$\Vert e_u  \Vert$'  +')-H1norm')
plt.title('H1 convergence test for velocity')
# plt.show()
plt.savefig(fuv_h1 + 'convergence_' + sur_fix + '.pdf', bbox_inches='tight')

plt.subplots(figsize=(11.5, 7))
for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap):
    temp_y = [np.log(eta_l2[i, i_step]) for i in range(npoint)]
    plt.plot(temp_x, temp_y, "*-", markersize = 12, linewidth = 4, color = 'mediumseagreen')
plt.legend(['slope=' + str(round(average_slope(eta_l2[:, i_step], hsize), 4)) for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap)], loc='upper left')
plt.grid(which='both')
plt.xlim([-4, 1])
plt.xticks(np.linspace(-4, 1, 6))
plt.ylim([-11, -1])
plt.yticks(np.linspace(-11, -1, 6))
plt.xlabel('log(h)')
plt.ylabel('log('+r'$\Vert e_{\eta} \Vert$'+')-L2norm')
plt.title('L2 convergence test for surface elevation')
# plt.show()
plt.savefig(feta_l2 + 'convergence_' + sur_fix + '.pdf', bbox_inches='tight')

plt.subplots(figsize=(11.5, 7))
for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap):
    temp_y = [np.log(eta_h1[i, i_step]) for i in range(npoint)]
    plt.plot(temp_x, temp_y, "*-", markersize = 12, linewidth = 4, color = 'mediumseagreen')
plt.legend(['slope=' + str(round(average_slope(eta_h1[:, i_step], hsize), 4)) for i_step in range(time_step_to_plot_start, time_step_to_plot_end, time_step_to_plot_gap)], loc='upper left')
plt.grid(which='both')
plt.xlim([-4, 1])
plt.xticks(np.linspace(-4, 1, 6))
plt.ylim([-11, -1])
plt.yticks(np.linspace(-11, -1, 6))
plt.xlabel('log(h)')
plt.ylabel('log('+r'$\Vert e_{\eta} \Vert$'+')-H1norm')
plt.title('H1 convergence test for surface elevation')
# plt.show()
plt.savefig(feta_h1 + 'convergence_' + sur_fix + '.pdf', bbox_inches='tight')
