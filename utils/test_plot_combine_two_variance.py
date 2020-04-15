import numpy as np
import sys
sys.path.append('/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src/')
from fenics import *


name = "test4"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
time_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
mesh_file = "inlet_adh_sswm_finer.xml"
#mesh_file = "Gulf_wind.xml"
time_file = "time_stamp_at_every_time_step.npy"
variance_u_input = "u_used_for_read_back_variance_in_domain.h5"
variance_v_input = "v_used_for_read_back_variance_in_domain.h5"
variance_combined = "uv_combined_variance_in_domain.xdmf"
variance_combined_h5file = "uv_combined_used_for_read_back_variance_in_domain.h5"

time_step = 500

mesh = Mesh(mesh_dir + mesh_file)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
D = FunctionSpace(mesh, "CG", 2)
var_u_f = Function(D)
var_v_f = Function(D)
var_uv_f = Function(D)

u_input = HDF5File(mesh.mpi_comm(), input_dir + variance_u_input, "r")
v_input = HDF5File(mesh.mpi_comm(), input_dir + variance_v_input, "r")
variance_uv_output = XDMFFile(output_dir + variance_combined)
variance_uvh5_output = HDF5File(mesh.mpi_comm(), output_dir + variance_combined_h5file, "w")


t = np.load(time_dir + time_file)

for k in range(time_step):
    dataset_u = "UVelocity/vector_%d"%k
    dataset_v = "VVelocity/vector_%d"%k
    u_input.read(var_u_f, dataset_u)
    v_input.read(var_v_f, dataset_v)
    var_uv_f.vector()[:] = var_u_f.vector()[:] + var_v_f.vector()[:]
    var_uv_f.rename("uv_var","variance u + variance v")

    if np.isnan(var_u_f.vector().array()).any():
        import pdb; pdb.set_trace()
    
    variance_uv_output.write(var_uv_f, float(t[k]))
    variance_uvh5_output.write(var_uv_f, "UandVvariance", float(t[k]))
    

    print "time step:", str(k), "done."

