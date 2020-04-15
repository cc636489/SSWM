

import dolfin

mesh = dolfin.UnitSquareMesh(20,20)
Q = dolfin.FunctionSpace(mesh, "CG", 2, 2)
F = dolfin.Function(Q)

hdf = dolfin.HDF5File(mesh.mpi_comm(), "a.h5", "w")

for i in range(10):
    hdf.write(F, "fun", float(dolfin.Constant(i)))

# Delete HDF5File object, closing file
del hdf


# Read functions back in
hdf = dolfin.HDF5File(mesh.mpi_comm(), "a.h5", "r")
attr = hdf.attributes("fun")
nsteps = attr['count']
for i in range(nsteps):
    dataset = "fun/vector_%d"%i
    attr = hdf.attributes(dataset)
    print 'Retrieving time step:', attr['timestamp']
    hdf.read(F, dataset)
