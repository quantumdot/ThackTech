import sys
import unittest
from ThackTech.Pipelines.ModuleParameter import ModuleParameter


def suite():
    return unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])


class TestModuleParameter(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_value_coercion(self):
        p = ModuleParameter('test', int, 0, value='3')
        self.assertEqual(p.value, 3)
        


if __name__ == "__main__":
    unittest.main()