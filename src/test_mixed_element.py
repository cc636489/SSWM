from fenics import *
import numpy as np

mesh = UnitSquareMesh(32,32)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
def boundaries(x,on_boundary):
    return on_boundary

bcs = DirichletBC(V, Constant(1.0), boundaries)

a = inner(grad(u), grad(v))*dx
L = Constant(1.0)*v*dx

u0 = Function(V)
solve(a==L, u0, bcs)

plot(u0)
interactive()

