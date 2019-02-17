from fenics import *
import numpy as np
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test3_res_wiLES/'
commonfolder = '20181001_generalized/'
modestatus = 'convergence_test/Initial_time_step_0.01_second/'
etafile_prefix = 'error_eta_timestep1_woles'
uvfile_prefix = 'error_u_timestep1_woles'

npoint = 4

feta = inputdir+testcase+commonfolder+modestatus+etafile_prefix
fuv = inputdir+testcase+commonfolder+modestatus+uvfile_prefix
eta, u = [], []

with open(feta, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        eta.append(float(line[0]))

with open(fuv, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        u.append(float(line[0]))

# hsize
hsize = [1./4., 1./8., 1./16., 1./32.]
tempx = [np.log(hsize[i]) for i in range(npoint)]


def average_slope(field, hsize):
    l = len(field)
    sum = 0.0
    for i in range(l-1):
        sum += (np.log(field[i]) - np.log(field[i+1])) / (np.log(hsize[i])-np.log(hsize[i+1]))
    return sum /(l-1)

# plots
surfix = 'initial_time'
plt.figure()
tempy = [np.log(u[i]) for i in range(npoint)]
plt.plot(tempx, tempy, "*-")
print "u:", average_slope(u, hsize)
plt.legend(['p='+str(round(average_slope(u, hsize), 4))], loc='best')
plt.grid(which='both')
plt.xlabel('log(h)')
plt.ylabel('log(|error_u|)')
# plt.xlim([-4, -0.5])
# plt.ylim([-4, -0.5])
plt.show()
# plt.savefig(fuv+'convergence_'+surfix+'_woles.png')

plt.figure()
tempy = [np.log(eta[i]) for i in range(npoint)]
plt.plot(tempx, tempy, "*-")
print "eta:", average_slope(eta, hsize)
plt.legend(['p='+str(round(average_slope(eta, hsize), 4))], loc='best')
plt.grid(which='both')
plt.xlabel('log(h)')
plt.ylabel('log(|error_eta|)')
# plt.xlim([-5, -0.5])
# plt.ylim([-5, -0.5])
plt.show()
# plt.savefig(feta+'convergence_'+surfix+'_woles.png')
