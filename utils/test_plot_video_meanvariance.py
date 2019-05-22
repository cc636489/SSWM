import numpy as np
import sys
sys.path.append('/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src/')
from fenics import *
from make_sto_basis import make_sto_basis

name = "IKE"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_stochastic/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_results/"
time_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"+name+"_bins/"
mesh_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
#mesh_file = "inlet_adh_sswm_finer.xml"
#u_file = "u_used_for_read_back_" + name + "_stochastic_0.2_0.75_0.13les_"
#eta_file = "eta_used_for_read_back_" + name + "_stochastic_0.2_0.75_0.13les_"
mesh_file = "Gulf_wind.xml"
u_file = "u_used_for_read_back_gulf_winds_ike_stochastic_"
eta_file = "eta_used_for_read_back_gulf_winds_ike_stochastic_"
time_file = "time_stamp_at_every_time_step.npy"
variance_eta_file = "eta_variance_in_domain.xdmf"
variance_u_file = "u_variance_in_domain.xdmf"
variance_v_file = "v_variance_in_domain.xdmf"
variance_eta_h5file = "eta_used_for_read_back_variance_in_domain.h5"
variance_u_h5file = "u_used_for_read_back_variance_in_domain.h5"
variance_v_h5file = "v_used_for_read_back_variance_in_domain.h5"

dist_name = "uniform"
sto_poly_deg = 1
sto_poly_dim = 1
coefficient = [0.8, 1.2]

basis = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
orth = basis["basis"]
n_modes = basis["n_modes"]

time_step = 500

mesh = Mesh(mesh_dir + mesh_file)
#mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
#mesh = RectangleMesh(Point(0, 0), Point(100, 50), 20, 10)
B = FunctionSpace(mesh, "CG", 1)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
D = FunctionSpace(mesh, "CG", 2)
var_u_f = Function(D)
var_v_f = Function(D)
var_eta_f = Function(B)

eta, eta_f, u, u_f = [], [], [], []
for mode in range(1, n_modes):
    eta.append(Function(B))
    u.append(Function(C))
    eta_input_file = eta_file + '{:02d}'.format(mode) + ".h5"
    eta_f.append(HDF5File(mesh.mpi_comm(), input_dir + eta_input_file, "r"))
    u_input_file = u_file + '{:02d}'.format(mode) + ".h5"
    u_f.append(HDF5File(mesh.mpi_comm(), input_dir + u_input_file, "r"))

variance_u = XDMFFile(output_dir + variance_u_file)
variance_v = XDMFFile(output_dir + variance_v_file)
variance_eta = XDMFFile(output_dir + variance_eta_file)
variance_etah5 = HDF5File(mesh.mpi_comm(), output_dir + variance_eta_h5file, "w")
variance_uh5 = HDF5File(mesh.mpi_comm(), output_dir + variance_u_h5file, "w")
variance_vh5 = HDF5File(mesh.mpi_comm(), output_dir + variance_v_h5file, "w")

t = np.load(time_dir + time_file) 

for k in range(time_step):

    dataset_u = "WaterVelocity/vector_%d"%k
    dataset_eta = "SurfaceElevation/vector_%d"%k

    tmpu, tmpv, tmpe = 0, 0, 0
    for mode in range(n_modes - 1):
        u_f[mode].read(u[mode], dataset_u)
        eta_f[mode].read(eta[mode], dataset_eta)
        m, n = u[mode].split(deepcopy = True)
        tmpu += m.vector().array() ** 2
        tmpv += n.vector().array() ** 2
        tmpe += eta[mode].vector().array() ** 2

    assert(all([tmpu[i] >= 0 for i in range(len(tmpu))]))
    assert(all([tmpe[i] >= 0 for i in range(len(tmpe))]))
    assert(all([tmpv[i] >= 0 for i in range(len(tmpv))]))

    var_u_f.vector()[:] = tmpu[:] 
    var_v_f.vector()[:] = tmpv[:]
    var_eta_f.vector()[:] = tmpe[:]

    var_u_f.rename("u_var", "variance of u")
    var_v_f.rename("v_var", "variance of v")
    var_eta_f.rename("eta_var", "variance of eta")

    if np.isnan(var_u_f.vector().array()).any():
        import pdb; pdb.set_trace()

    variance_u.write(var_u_f, float(t[k]))
    variance_v.write(var_v_f, float(t[k]))
    variance_eta.write(var_eta_f, float(t[k]))

    variance_uh5.write(var_u_f, "UVelocity", float(t[k]))
    variance_vh5.write(var_v_f, "VVelocity", float(t[k]))
    variance_etah5.write(var_eta_f, "SurfaceElevation", float(t[k]))

    print "time step:", str(k), "done."


