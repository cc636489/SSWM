from dolfin import *

#parameters.form_compiler.representation = "quadrature"  # Works
#parameters.form_compiler.representation = "tensor"      # Works
#parameters.form_compiler.representation = "uflacs"      # Fails
#parameters.form_compiler.representation = "tsfc"        # Fails

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'P', 1)
v = TestFunction(V)

element = FiniteElement("CG", mesh.ufl_cell(), 1)
B = FunctionSpace(mesh, element)
p = Expression('2.71 * x[0]', element=B.ufl_element())
assemble(p.dx(0)*v*dx)
