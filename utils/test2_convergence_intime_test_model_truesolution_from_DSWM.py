" Write a nonlinear shallow water equation solver."


from fenics import *
from mshr import *
import numpy as np
import sympy as sp

truth_step = 256
loc = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test2_res/convergence_intime_test_result/"
#u_truth_filename = "u_used_for_read_back_test2_timestep_0.0048828125_truth_00.h5"
#eta_truth_filename = "eta_used_for_read_back_test2_timestep_0.0048828125_truth_00.h5"
u_truth_filename = "u_used_for_read_back_test2_timestep_0.0390625_truth_00.h5"
eta_truth_filename = "eta_used_for_read_back_test2_timestep_0.0390625_truth_00.h5"

x0 = 0
y0 = 0
x1 = 100
y1 = 50
nx = 200 
ny = 100
start_time = 0.
end_time = 10.

timestepsize = [5,2.5,1.25,0.625,0.3125]

# build up mesh
mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), nx, ny)

# Get Funciton space
V = VectorFunctionSpace(mesh, 'CG', 2, dim=2)
Q = FunctionSpace(mesh, 'CG', 1)
u_truth = Function(V, name="u_truth")
u_model = Function(V, name="u_model")
eta_truth = Function(Q, name="eta_truth")
eta_model = Function(Q, name="eta_model")

# extract true solution
u_f = HDF5File(mesh.mpi_comm(), loc + u_truth_filename, "r")
eta_f = HDF5File(mesh.mpi_comm(), loc + eta_truth_filename, "r")
dataset_u = "WaterVelocity/vector_%d"%truth_step
dataset_eta = "SurfaceElevation/vector_%d"%truth_step
u_f.read(u_truth, dataset_u)
eta_f.read(eta_truth, dataset_eta)

# calculate which timestep to extract
timestep_to_extract = [int((end_time - start_time)/x) for x in timestepsize]

# calculate the error norms.
for i, scale in enumerate([1,2,4,8,16]):
    u_file = HDF5File(mesh.mpi_comm(), loc + u_truth_filename[:36]+str(timestepsize[i])+u_truth_filename[-6:], "r")
    eta_file = HDF5File(mesh.mpi_comm(), loc + eta_truth_filename[:38]+str(timestepsize[i])+eta_truth_filename[-6:], "r")
    data_u = "WaterVelocity/vector_%d"%timestep_to_extract[i]
    data_eta = "SurfaceElevation/vector_%d"%timestep_to_extract[i]
    u_file.read(u_model, data_u)
    eta_file.read(eta_model, data_eta)
    error_u_l2 = errornorm(u_truth, u_model, "L2", 8)
    error_eta_l2 = errornorm(eta_truth, eta_model, "L2", 4)
    f = open(loc+'error_u_l2_'+str(scale), "w+")
    f.write(str(error_u_l2))
    f.write("\n")
    f.close()
    f = open(loc+'error_eta_l2_'+str(scale), "w+")
    f.write(str(error_eta_l2))
    f.write("\n")
    f.close()

