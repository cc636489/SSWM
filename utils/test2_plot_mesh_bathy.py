from fenics import *

dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test2_res/20181001_generalized/bathy_mesh_plot/"

mesh = RectangleMesh(Point(0, 0), Point(100, 50), 100, 50)
bathy = Constant(20)
V = FunctionSpace(mesh, 'CG', 1)
bathy = interpolate(bathy, V)
f = XDMFFile(dir+'test2_mesh_bathy.xdmf')
f.write(bathy)

