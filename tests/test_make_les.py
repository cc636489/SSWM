

import unittest
from src.make_les import LES
from src.make_mesh import make_mesh
from fenics import VectorFunctionSpace, FunctionSpace, Constant, interpolate


class MakeLesTestCase(unittest.TestCase):

    def setUp(self):
        domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
        mesh = make_mesh(domain)
        self.v = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        self.b = FunctionSpace(mesh, "CG", 1)
        exp = Constant([10.0, 0.0])
        self.u = interpolate(exp, self.v)

    def test_les(self):
        tmp = LES(self.b, self.u, 1.0)
        tmp.solve()
        self.assertTrue(1)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
