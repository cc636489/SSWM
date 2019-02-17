from fenics import *
import numpy as np
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test2_res/'
commonfolder = '20181001_generalized/'
modestatus = 'mode_0_analytical_comparison/'
etafile = 'eta_sloth_25_25.csv'
uvfile = 'uv_sloth_25_25.csv'

feta = inputdir+testcase+commonfolder+modestatus+etafile
fuv = inputdir+testcase+commonfolder+modestatus+uvfile
time, eta, u, v = [], [], [], []

with open(feta, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
            print len(line.split(','))
            print line.split(',')
            continue
        line = line.split(',')
        eta.append(float(line[0]))
        time.append(float(line[1]))

with open(fuv, 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
            print len(line.split(','))
            print line.split(',')
            continue
        line = line.split(',')
        u.append(float(line[0]))
        v.append(float(line[1]))

# analytical solution.
x = 25.0
L = 100.0
H = 20.0
a = 0.1
g = 9.81
pi = 3.14159

eta_true = list(map(lambda t: a*cos(pi*x/L)*cos(pi*sqrt(g*H)/L*t), time))
u_true = list(map(lambda t: a*sqrt(g*H)/H*sin(pi*x/L)*sin(pi*sqrt(g*H)/L*t), time))
v_true = list(map(lambda _: 0, time))

# plots
print time
print [eta[i]-eta_true[i] for i in range(len(eta))]
print [u[i]-u_true[i] for i in range(len(u))]
print [v[i]-v_true[i] for i in range(len(v))]

plt.figure(1)
plt.plot(time, eta, '*-', time, eta_true, '.-')
plt.legend(['model solution', 'analytical solution'], loc='best')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(sec)')
plt.ylabel('eta(m)')
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
plt.grid()
# plt.show()
plt.savefig(inputdir+testcase+commonfolder+modestatus+'eta_line_comparison_25_25.png')

plt.figure(2)
plt.plot(time, u, '*-', time, u_true, '.-')
plt.legend(['model solution', 'analytical solution'], loc='best')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(seconds)')
plt.ylabel('u(m/s)')
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
plt.grid()
# plt.show()
plt.savefig(inputdir+testcase+commonfolder+modestatus+'u_line_comparison_25_25.png')


plt.figure(3)
plt.plot(time, v, '*-', time, v_true, '.-')
plt.legend(['model solution', 'analytical solution'], loc='best')
plt.ylim([min(eta)*1.5, max(eta)*1.5])
plt.xlabel('time(seconds)')
plt.ylabel('v(m/s)')
plt.grid()
plt.ylim([-0.1,0.1])
plt.xlim([0,50])
# plt.show()
plt.savefig(inputdir+testcase+commonfolder+modestatus+'v_line_comparison_25_25.png')


