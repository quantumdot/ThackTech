import unittest


def suite():
    suites = []
    
    import test_FileInfo
    suites.append(test_FileInfo.suite())
    
    import test_Context
    suites.append(test_Context.suite())
    
    import test_ModuleParameter
    suites.append(test_ModuleParameter.suite())

    alltests = unittest.TestSuite(suites)
    return alltests 
