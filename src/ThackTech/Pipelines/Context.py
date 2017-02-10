



class ModuleContext(object):

    def __init__(self, pipeline_name, step_num, module_name):
        self.__pipeline = pipeline_name
        self.__step = step_num
        self.__module = module_name
    
    @property
    def pipeline(self):
        return self.__pipeline
    
    @property
    def step(self):
        return self.__step
    
    @property
    def module(self):
        return self.__module
    
#end class ModuleContext

class ModuleRunContext(ModuleContext):

    def __init__(self, pipeline_name, step_num, module_name, logfile, sample):
        super(ModuleContext, self).__init__(pipeline_name, step_num, module_name)
        self.__logfile = logfile
        self.__sample = sample
        
    @property
    def sample(self):
        return self.__sample
    
    @property
    def logfile(self):
        return self.__logfile
    
#end class ModuleRunContext

