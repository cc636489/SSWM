

from fenics import *
import numpy as np
import subprocess
import sys
sys.path.insert(0, '/Users/chenchen/gloria/test_fenics/SupportingRun/test5_SpurDike/')
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
      'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick
import pandas as pd

name = "test5"
#truthname = "Base"
input_dir = "/Users/chenchen/gloria/test_fenics/SupportingRun/test5_SpurDike/"
output_dir = "/Users/chenchen/gloria/test_fenics/SupportingRun/test5_SpurDike/"
mesh_dir = "/Users/chenchen/gloria/test_fenics/SupportingRun/test5_SpurDike/"
mesh_file = "dike.xml"
u_file = "u_used_for_read_back_"+name.lower()+"_generalized_"
eta_file = "eta_used_for_read_back_"+name.lower()+"_generalized_"
#truth_file = "Near"+truthname+"_"
number = "1.5_0.0015"
umin = -0.2
umax = 0.4
form = ".png"

m = 14
n = 0.152
test_node_x = [m+2*n, m+4*n, m+6*n, m+8*n]
test_node_y = np.linspace(0, 0.92, 50)
truthselect = [2, 4, 6, 8]


k = 160

mesh = Mesh(mesh_dir + mesh_file)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)

eta = Function(B)
u = Function(C)
eta_input_file = eta_file + number +"00.h5"
u_input_file = u_file + number + "00.h5"
eta_f = HDF5File(mesh.mpi_comm(), input_dir + eta_input_file, "r")
u_f = HDF5File(mesh.mpi_comm(), input_dir + u_input_file, "r")

print("start: timestep = " + str(k))

dataset_u = "WaterVelocity/vector_%d"%k
dataset_eta = "SurfaceElevation/vector_%d"%k

u_f.read(u, dataset_u)
eta_f.read(eta, dataset_eta)


for i, x in enumerate(test_node_x): 

    print("    test_node_x: x = " + str(x)) 

    u_list = [u(x, y)[0] for y in test_node_y]
    
    matrix = pd.read_csv(input_dir+"NearSurface_"+str(truthselect[i])+".csv", sep=',', header=None)

    u_truth_surface = matrix.iloc[:, 0].tolist()
    y_truth_surface = matrix.iloc[:, 1].tolist()

    matrix = pd.read_csv(input_dir+"NearBase_"+str(truthselect[i])+".csv", sep=',', header=None)

    u_truth_Base = matrix.iloc[:, 0].tolist()
    y_truth_Base = matrix.iloc[:, 1].tolist()

    #import pdb;pdb.set_trace()
    plt.figure(figsize=[7, 7])
    plt.ylim(0, 0.7)
    plt.ylabel('y(m)')
    plt.xlim(umin,umax)
    plt.xticks(np.linspace(umin, umax, 4))
    plt.xlabel('x-direction velocity(m/s)')
    #plt.title('x-direction water velocity profile ' + str(k * step_size) + ' sec')
    plt.plot(u_list, test_node_y, '*-', color='orange')
    plt.scatter(u_truth_surface, y_truth_surface, marker='o', s=60,  facecolors='none', edgecolors='r')
    plt.scatter(u_truth_Base, y_truth_Base, marker='o', s=60,  facecolors='none', edgecolors='b')
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
    plt.legend(['DSWM','z=0.85H Data','z=0.03H Data'])
    plt.savefig(output_dir + 'u_all_profile_' + number + str(i) +'_timestep_'+str(k)+form)
    plt.close()
