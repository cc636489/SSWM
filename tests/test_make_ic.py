

import unittest
from src.make_mesh import make_mesh
from src.make_sto_basis import make_sto_basis
from src.make_ic import make_initial_object_list
from fenics import FiniteElement, VectorElement, FunctionSpace, MixedElement


class MakeStoTestCase(unittest.TestCase):

    def setUp(self):
        mesh = make_mesh({"rectangle": (0., 0., 100., 50., 20, 10)})
        self.basis_str = make_sto_basis("uniform", 2, 2, [0.18, 0.22, -0.001, 0.001])
        n_modes = self.basis_str.get("n_modes")
        v = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
        q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
        string_v = "["
        string_q = "["
        for i in range(n_modes - 1):
            string_v = string_v + "v,"
            string_q = string_q + "q,"
        string_v = string_v + "v]"
        string_q = string_q + "q]"
        self.w = FunctionSpace(mesh, MixedElement(eval(string_v)))
        self.p = FunctionSpace(mesh, MixedElement(eval(string_q)))

    def test_none_ic_1(self):
        initial_u = {"flat": (0., 0.)}
        initial_eta = {"flat": 0.0}
        u_ic_function, eta_ic_function = make_initial_object_list(initial_u, initial_eta,
                                                                  self.basis_str, self.w, self.p, False)
        self.assertIsInstance(u_ic_function, object)
        self.assertIsInstance(eta_ic_function, object)

    def test_none_ic_2(self):
        initial_u = {"flat": (0., 0.)}
        initial_eta = {"vary": '20.*x*q1'}
        u_ic_function, eta_ic_function = make_initial_object_list(initial_u, initial_eta,
                                                                  self.basis_str, self.w, self.p, False)
        self.assertIsInstance(u_ic_function, object)
        self.assertIsInstance(eta_ic_function, object)

    def test_none_ic_3(self):
        initial_u = {"vary": ('2*sp.sin(x+y)*q0', '3*x**4*q1')}
        initial_eta = {"flat": 0.0}
        u_ic_function, eta_ic_function = make_initial_object_list(initial_u, initial_eta,
                                                                  self.basis_str, self.w, self.p, False)
        self.assertIsInstance(u_ic_function, object)
        self.assertIsInstance(eta_ic_function, object)

    def test_none_ic_4(self):
        initial_u = {"vary": ('2*sp.sin(x+y)*q0', '3*x**4*q1')}
        initial_eta = {"vary": '20.*x*y*q0*q1'}
        u_ic_function, eta_ic_function = make_initial_object_list(initial_u, initial_eta,
                                                                  self.basis_str, self.w, self.p, False)
        self.assertIsInstance(u_ic_function, object)
        self.assertIsInstance(eta_ic_function, object)

    def test_none_ic_5(self):
        initial_u = {"flat": (0.0, 0.0)}
        initial_eta = {"flat": 0.0}
        u_ic_function, eta_ic_function = make_initial_object_list(initial_u, initial_eta,
                                                                  self.basis_str, self.w, self.p, False)
        self.assertIsInstance(u_ic_function, object)
        self.assertIsInstance(eta_ic_function, object)
