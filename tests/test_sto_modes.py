

import unittest
from src.make_sto_basis import make_sto_basis
from src.make_sto_modes import make_sto_modes


class MakeStoModesTestCase(unittest.TestCase):

    def setUp(self):
        dist_name = "uniform"
        sto_poly_deg = 2
        sto_poly_dim = 2
        coefficient = [0.18, 0.22, -0.001, 0.001]
        self.basis_str = make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, coefficient)

    def test_1(self):
        tmp = make_sto_modes(self.basis_str, "0.0")
        self.assertEqual(tmp, ['0', '0', '0', '0', '0', '0'], "edge case 0.0, generation error")

    def test_2_1(self):
        tmp = make_sto_modes(self.basis_str, "1.0")
        self.assertAlmostEqual(float(tmp[0]), 1.0, places=7, msg="edge case 1.0, generation error")

    def test_2_2(self):
        tmp = make_sto_modes(self.basis_str, "1.0")
        n_modes = self.basis_str.get("n_modes")
        self.assertAlmostEqual(float(tmp[n_modes-1]), 0.0, places=7, msg="edge case 1.0 generation error")

    def test_3(self):
        tmp = make_sto_modes(self.basis_str, "2*x*q0+3*y*q1")
        n_modes = self.basis_str.get("n_modes")
        self.assertIsInstance(tmp[n_modes-1], str, msg="string case generation error")

    def test_4(self):
        tmp = make_sto_modes(self.basis_str, "2*x*q0+3*y*q1", "q0*q1")
        self.assertIsInstance(tmp[0], list, msg="string case generation error")
        self.assertIsInstance(tmp[1], list, msg="string case generation error")

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
