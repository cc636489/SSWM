

import unittest
from src.make_sto_ijk import make_sto_ijk
from numpy import ndarray
import chaospy as cp
from src.make_sto_basis import make_sto_basis


class MakeTestCase(unittest.TestCase):

    def setUp(self):
        self.dist_name = "uniform"
        self.sto_poly_deg = 2
        self.sto_poly_dim = 2
        self.coefficient = [0.18, 0.22, -0.001, 0.001]

    def test_make_sto_basis(self):
        basis_str = make_sto_basis(self.dist_name, self.sto_poly_deg, self.sto_poly_dim, self.coefficient)
        ort = basis_str.get('basis')
        n_modes = basis_str.get('n_modes')
        cdf = basis_str.get('joint_cdf')
        for i in range(n_modes):
            for j in range(n_modes):
                if i == j:
                    self.assertTrue(abs(cp.E(ort[i] * ort[j], cdf) - 1.0) <= 1.e-10,
                                    "Not satisfy orthogonality condition for basis " + str(i) + " and basis " + str(j))
                else:
                    self.assertTrue(abs(cp.E(ort[i] * ort[j], cdf)) <= 1.e-10,
                                    "Not satisfy orthogonality condition for basis " + str(i) + " and basis " + str(j))
