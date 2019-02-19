

import unittest
from src.make_sto_basis import make_sto_basis
from src.make_mesh import make_mesh
from src.make_bc import make_boundary_object_list
from fenics import FiniteElement, VectorElement, FunctionSpace, MixedElement


class MakeStoTestCase(unittest.TestCase):

    def setUp(self):
        self.domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
        self.mesh = make_mesh(self.domain)
        self.basis_str = make_sto_basis("uniform", 2, 2, [0.18, 0.22, -0.001, 0.001])
        n_modes = self.basis_str.get("n_modes")
        v = VectorElement(family='CG', cell=self.mesh.ufl_cell(), degree=2, dim=2)
        q = FiniteElement(family='CG', cell=self.mesh.ufl_cell(), degree=1)
        string_v = "["
        string_q = "["
        for i in range(n_modes - 1):
            string_v = string_v + "v,"
            string_q = string_q + "q,"
        string_v = string_v + "v]"
        string_q = string_q + "q]"
        self.w = FunctionSpace(self.mesh, MixedElement(eval(string_v)))
        self.p = FunctionSpace(self.mesh, MixedElement(eval(string_q)))

    def test_none_bc_1(self):
        bc_file = None
        boundary_u = {1: ("sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0", "q1"), 3: "free_slip_in_x", 4: "no_slip"}
        boundary_eta = {2: 0.0}
        u_bc_object_list, eta_bc_object_list, u_time_dependent, \
            eta_time_dependent, u_list_expression, eta_list_expression = \
            make_boundary_object_list(boundary_u, boundary_eta, self.basis_str, bc_file,
                                      self.mesh, self.domain, self.w, self.p, 0, 1, 1)
        self.assertIsInstance(u_bc_object_list, list)
        self.assertIsInstance(eta_bc_object_list, list)
        self.assertIsInstance(u_list_expression, object)
        self.assertIsInstance(eta_list_expression, object)
        self.assertTrue(u_time_dependent)
        self.assertFalse(eta_time_dependent)

    def test_none_bc_2(self):
        bc_file = None
        boundary_u = {1: (10., 0.), 3: "no_slip", 4: "free_slip_in_x"}
        boundary_eta = {2: "M2 special"}
        u_bc_object_list, eta_bc_object_list, u_time_dependent, \
            eta_time_dependent, u_list_expression, eta_list_expression = \
            make_boundary_object_list(boundary_u, boundary_eta, self.basis_str, bc_file,
                                      self.mesh, self.domain, self.w, self.p, 0, 1, 1)
        self.assertIsInstance(u_bc_object_list, list)
        self.assertIsInstance(eta_bc_object_list, list)
        self.assertIsInstance(u_list_expression, object)
        self.assertIsInstance(eta_list_expression, object)
        self.assertFalse(u_time_dependent)
        self.assertTrue(eta_time_dependent)

    def test_none_bc_3(self):
        bc_file = None
        boundary_u = {1: "free_slip_in_y", 3: "free_slip_in_x", 4: "no_slip"}
        boundary_eta = {2: "sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0"}
        u_bc_object_list, eta_bc_object_list, u_time_dependent, \
            eta_time_dependent, u_list_expression, eta_list_expression = \
            make_boundary_object_list(boundary_u, boundary_eta, self.basis_str, bc_file,
                                      self.mesh, self.domain, self.w, self.p, 0, 1, 1)
        self.assertIsInstance(u_bc_object_list, list)
        self.assertIsInstance(eta_bc_object_list, list)
        self.assertIsInstance(u_list_expression, object)
        self.assertIsInstance(eta_list_expression, object)
        self.assertFalse(u_time_dependent)
        self.assertTrue(eta_time_dependent)
