import sys
import unittest
from ThackTech.Pipelines.FileInfo import FileContext
from ThackTech.Pipelines.Context import BaseModuleContext, ModuleRunContext


def suite():
    return unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])


class TestContext(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        
    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def test_base_context_equality(self):
        c1 = BaseModuleContext('pipeline1', 1, 'module1')
        c2 = BaseModuleContext('pipeline1', 1, 'module1')
        self.assertTrue(c1 == c2)
        self.assertFalse(c1 != c2)
    
    def test_base_context_non_equality(self):
        c1 = BaseModuleContext('pipeline1', 1, 'module1')
        c2 = BaseModuleContext('pipeline1', 2, 'module2')
        self.assertFalse(c1 == c2)
        self.assertTrue(c1 != c2)
        
    def test_file_context_equality(self):
        c1 = FileContext('pipeline1', 1, 'module1', 'file1')
        c2 = FileContext('pipeline1', 1, 'module1', 'file1')
        self.assertTrue(c1 == c2)
        self.assertFalse(c1 != c2)
        
    def test_file_context_non_equality(self):
        c1 = FileContext('pipeline1', 1, 'module1', 'file1')
        c2 = FileContext('pipeline1', 1, 'module1', 'file2')
        self.assertFalse(c1 == c2)
        self.assertTrue(c1 != c2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
