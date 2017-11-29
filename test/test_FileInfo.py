import os
import sys
import unittest
from ThackTech.Pipelines.FileInfo import FileInfo


def suite():
    return unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])


class TestFileInfo(unittest.TestCase):

    def setUp(self):
        #self.ntf = NamedTemporaryFile()
        unittest.TestCase.setUp(self)
        
    def tearDown(self):
        #self.ntf.close()
        unittest.TestCase.tearDown(self)

    def test_fileinfo_properties(self):
        finfo = FileInfo('/path/to/some/file.ext1.ext2')
        self.assertEqual(finfo.basename, 'file.ext1.ext2')
        self.assertEqual(finfo.basename_with_ext('ext3'), 'file.ext1.ext3')
        self.assertEqual(finfo.dirname, os.path.abspath('/path/to/some'))
        self.assertEqual(finfo.ext, '.ext2')
        self.assertEqual(finfo.filename, 'file.ext1')
        self.assertEqual(finfo.filename_strip_all_ext, 'file')
        

    def test_fileinfo_attributes(self):
        attributes = {'attr1': 1, 'attr2': 2}
        finfo = FileInfo('/path/to/some/file.ext1.ext2', **attributes)
        self.assertDictEqual(finfo.attributes, attributes)
        self.assertTrue(finfo.has_attribute('attr1'))
        self.assertTrue(finfo.has_attribute_value('attr1', 1))
        self.assertFalse(finfo.has_attribute_value('attr1', 5))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()