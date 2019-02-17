from fenics import *

mesh = RectangleMesh(Point(0, 0), Point(1000, 200), 100, 20)
Q = FunctionSpace(mesh, 'CG', 1)


class bottomExpression(Expression):

    def eval(self, value, x):
        if x[0] < 400 or x[0] > 600:
            value[0] = 5.0
        else:
            value[0] = -3*((x[0]-500)/100)**4 + 6*((x[0]-500)/100)**2 + 2

    def value_shape(self):
        return (1, )


bathy = bottomExpression(element=Q.ufl_element())
bathy_f = interpolate(bathy, Q)

f = XDMFFile('test3_mesh_bathy.xdmf')
f.write(bathy_f)