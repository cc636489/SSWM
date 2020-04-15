from fenics import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal', 'size'   : 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test2_res/'
modestatus = 'ETA_U_SIMULATION_RESULTS/'
point_x = 75
point_y = 25
etafile = 'eta_sloth_' + str(point_x) + '_' + str(point_y) + '.csv'
uvfile = 'uv_sloth_' + str(point_x) + '_' + str(point_y) + '.csv'

feta = inputdir+testcase+modestatus+etafile
fuv = inputdir+testcase+modestatus+uvfile
time, eta, u, v = [], [], [], []

with open(feta, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
#            print len(line.split(','))
#            print line.split(',')
            continue
        line = line.split(',')
        eta.append(float(line[0]))
        time.append(float(line[2]))

with open(fuv, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
#            print len(line.split(','))
#            print line.split(',')
            continue
        line = line.split(',')
        u.append(float(line[0]))
        v.append(float(line[1]))

# analytical solution.
x = point_x 
L = 100.0
H = 20.0
a = 0.1
g = 9.81
pi = 3.14159

eta_true = list(map(lambda t: a*cos(pi*x/L)*cos(pi*sqrt(g*H)/L*t), time))
u_true = list(map(lambda t: a*sqrt(g*H)/H*sin(pi*x/L)*sin(pi*sqrt(g*H)/L*t), time))
v_true = list(map(lambda _: 0, time))

# plots
#print time
#print [eta[i]-eta_true[i] for i in range(len(eta))]
#print [u[i]-u_true[i] for i in range(len(u))]
#print [v[i]-v_true[i] for i in range(len(v))]

plt.subplots(figsize=(13,7))
plt.plot(time, eta,'-', color = 'mediumseagreen', linewidth = 4)
plt.plot(time, eta_true,'*', color= 'darkorange', markersize = 12)
plt.legend(['DSWM', 'Analytical'], loc='upper right')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(sec)')
plt.ylabel('eta(m)')
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
plt.grid()
# plt.show()
plt.savefig(inputdir+testcase+modestatus+'eta_line_comparison_' + str(point_x) + '_' + str(point_y) + '.pdf', bbox_inches='tight')

plt.subplots(figsize=(13,7))
plt.plot(time, u,'-', color = 'mediumseagreen', linewidth = 4)
plt.plot(time, u_true,'*', color= 'darkorange', markersize = 12)
plt.legend(['DSWM', 'Analytical'], loc='upper right')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(seconds)')
plt.ylabel('u(m/s)')
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
plt.grid()
# plt.show()
plt.savefig(inputdir+testcase+modestatus+'u_line_comparison_' + str(point_x) + '_' + str(point_y) + '.pdf', bbox_inches='tight')


plt.subplots(figsize=(13,7))
plt.plot(time, v,'-', color = 'mediumseagreen', linewidth = 4)
plt.plot(time, v_true,'*', color= 'darkorange', markersize = 12)
plt.legend(['DSWM', 'Analytical'], loc='upper right')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(seconds)')
plt.ylabel('v(m/s)')
plt.grid()
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
# plt.show()
plt.savefig(inputdir+testcase+modestatus+'v_line_comparison_' + str(point_x) + '_' + str(point_y) + '.pdf',bbox_inches='tight')


