

import unittest
from src.make_mesh import make_mesh


class MakeMeshTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_builtin_mesh(self):
        domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
        tmp = make_mesh(domain)
        self.assertIsInstance(tmp, object, "builtin mesh generation error")

    def test_import_mesh(self):
        domain = {"importfile": "../input/Gulf_wind.xml"}
        tmp = make_mesh(domain)
        self.assertIsInstance(tmp, object, "generalized mesh generation error")

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
