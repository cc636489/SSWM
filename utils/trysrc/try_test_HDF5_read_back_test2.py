
from fenics import *

mesh = RectangleMesh(Point(0,0), Point(100, 50), 100, 50)
a = FunctionSpace(mesh, "CG", 1)
eta = Function(a)
f = HDF5File(mesh.mpi_comm(), "/workspace/Documentation/Research_Doc/"
                              "SFEM_Doc/4-NS-results-and-tests/regression_test/"
                              "eta_used_for_read_back_test2_generalized_00.h5", "r")
attr = f.attributes("SurfaceElevation")
n_steps = attr["count"]
for i in range(n_steps):
    dataset = "SurfaceElevation/vector_%d"%i
    attr = f.attributes(dataset)
    f.read(eta, dataset)
    plot(eta)
