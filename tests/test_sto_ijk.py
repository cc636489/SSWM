

import unittest
from src.make_sto_basis import make_sto_basis
from src.make_sto_ijk import make_sto_ijk


class MakeStoIJKTestCase(unittest.TestCase):

    def setUp(self):
        self.dist_name = "uniform"
        self.coefficient = [0.18, 0.22]

    def test_unity(self):
        basis_str = make_sto_basis(dist_name=self.dist_name, sto_poly_deg=0, sto_poly_dim=1, co_eff=self.coefficient)
        tmp = make_sto_ijk(basis_str)
        self.assertTrue(tmp[0][0][0] - 1.0 < 1e-10, "stoIJK calculation error")

    def test_orthogonality(self):
        basis_str = make_sto_basis(dist_name=self.dist_name, sto_poly_deg=2, sto_poly_dim=1, co_eff=self.coefficient)
        tmp = make_sto_ijk(basis_str)
        self.assertTrue(tmp[0][1][2] < 1e-10, "stoIJK calculation error")

    def tearDown(self):
        pass
