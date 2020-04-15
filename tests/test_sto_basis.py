

import unittest
import chaospy as cp
from src.make_sto_basis import make_sto_basis


class MakeStoBasisTestCase(unittest.TestCase):

    def setUp(self):
        dist_name = "uniform"
        sto_poly_deg = 2
        sto_poly_dim = 2
        coefficient = [0.18, 0.22, -0.001, 0.001]
        basis_str = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)
        self.ort = basis_str.get('basis')
        self.n_modes = basis_str.get('n_modes')
        self.cdf = basis_str.get('joint_cdf')

    def test_orthogonality(self):
        for i in range(self.n_modes):
            for j in range(i + 1, self.n_modes):
                self.assertTrue(abs(cp.E(self.ort[i] * self.ort[j], self.cdf)) <= 1.e-10,
                                "Not satisfy orthogonality condition for basis " + str(i) + " and basis " + str(j))

    def test_unity(self):
        for i in range(self.n_modes):
                self.assertTrue(abs(cp.E(self.ort[i] * self.ort[i], self.cdf) - 1.0) <= 1.e-10,
                                "Not satisfy unity condition for basis " + str(i))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
