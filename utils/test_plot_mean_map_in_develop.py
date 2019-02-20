import numpy as np
import sys
sys.path.append('/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-bitbucket/src/')
from makemesh import MakeMesh
from fenics import *
from dolfin import *



loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
dir = 'test2_res/'
name = 'eta'
file = loc + dir + name + '_fem_verified_test2_00.h5'


domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
# domain = {"rectangle": (0., 0., 1000., 200., 100, 20)}
mesh = MakeMesh(domain)
V = VectorFunctionSpace(mesh, 'CG', 2, 2)
Q = FunctionSpace(mesh, 'CG', 1)

u = Function(V)
eta = Function(Q)

f = XDMFFile('eta.xdmf')
f.write(eta)
f.close()

input_file = XDMFFile(mpi_comm_world(), 'eta.xdmf', "r")
input_file.read_checkpoint(u, '')
input_file.close()

plot(u)
