



class BaseModuleContext(object):

    def __init__(self, pipeline_name, step_num, module_name):
        self.__pipeline = pipeline_name
        self.__step = int(step_num)
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
    
    def __str__(self):
        return "%s_%s_%s" % (self.pipeline, self.step, self.module)
#end class ModuleContext




class ModuleRunContext(BaseModuleContext):

    def __init__(self, pipeline_name, step_num, module_name, logfile, sample):
        super(ModuleRunContext, self).__init__(pipeline_name, step_num, module_name)
        self.__logfile = logfile
        self.__cxt.sample = sample
        
    @property
    def sample(self):
        return self.__cxt.sample
    
    @property
    def log(self):
        return self.__logfile
    
    def __str__(self):
        return "%s_%s" % (super(ModuleRunContext, self).__str__(), self.sample.name)
#end class ModuleRunContext

