

import unittest
from make_mesh import make_mesh
from make_bath import make_bath
from make_sto_basis import make_sto_basis
from fenics import FiniteElement, FunctionSpace, MixedElement


class MakeStoTestCase(unittest.TestCase):

    def setUp(self):
        domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
        mesh = make_mesh(domain)
        self.basis_str = make_sto_basis("uniform", 2, 2, [0.18, 0.22, 0.9, 1.1])
        n_modes = self.basis_str.get("n_modes")
        string_q = "["
        for i in range(n_modes - 1):
            string_q = string_q + "q,"
        string_q = string_q + "q]"
        if n_modes == 1:
            self.p = FunctionSpace(mesh, "CG", 1)
        else:
            q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
            self.p = FunctionSpace(mesh, MixedElement(eval(string_q)))

    def test_flat_bath(self):
        bath = {"flat": 20.}
        tmp = make_bath(bath, self.basis_str, self.p)
        self.assertIsInstance(tmp, object, "flat bathymetry generation error")

    def test_vary_bath(self):
        bath = {"vary": '2*x*q0+3*y*q1'}
        tmp = make_bath(bath, self.basis_str, self.p)
        self.assertIsInstance(tmp, object, "spacial-varied 1 bathymetry generation error")

    def test_type_1_bath(self):
        bath = {"class": ['type1', 400, 600, '5.0', '-3.0*((x-500.0)/100.0)**4*q0 + 6*((x-500.)/100.)**2*q1 + 2.']}
        tmp = make_bath(bath, self.basis_str, self.p)
        self.assertIsInstance(tmp, object, "spacial-varied type1 bathymetry generation error")

    def test_type_2_bath(self):
        bath = {"class": ['type2', 2150, '-14.0/2150.0*x + 19.0', '5.0']}
        tmp = make_bath(bath, self.basis_str, self.p)
        self.assertIsInstance(tmp, object, "spacial-varied type2 bathymetry generation error")

    def test_type_3_bath(self):
        bath = {"class": ['type3', "../input/inlet_adh_sswm_finer.nc"]}
        tmp = make_bath(bath, self.basis_str, self.p)
        self.assertIsInstance(tmp, object, "spacial-varied type3 bathymetry generation error")


if __name__ == '__main__':
    unittest.main()
