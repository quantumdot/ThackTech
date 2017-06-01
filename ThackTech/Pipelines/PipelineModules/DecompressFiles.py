from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class DecompressFiles(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='DecompressFiles', short_description='Decompress Files')
		super_args.update(**kwargs)
		super(DecompressFiles, self).__init__(**super_args)
	#end __init__()
	
	def __declare_parameters(self):
		self.add_parameter(ModuleParameter('overwrite',	bool, False, desc="Overwrite the decompressed file if it already exists."))
	#end __declare_parameters()
	
	def __declare_resolvers(self):
		self._name_resolver('files')
	#end __declare_resolvers()

	
	def run(self, cxt):
		if not cxt.sample.has_attribute('decompressed_files'):
			cxt.sample.set_attribute('decompressed_files', [])
		files_to_decompress = self.resolve_input('files', cxt)
		decompress_procs = []
		
		for f in files_to_decompress:
			if filetools.is_compressed(f.fullpath):
				cxt.log.write("\t-> Decompressing file {}....\n".format(f.fullpath))
				cxt.log.flush()
				d, p = filetools.extract(f.fullpath, cxt.sample.dest, overwrite=self.get_parameter_value('overwrite'), stderr=cxt.log)
				cxt.sample.set_attribute('decompressed_files', cxt.sample.get_attribute('decompressed_files') + [(d, f.fullpath)])
				f._set_path(d)
				decompress_procs.append(p)
			else:
				cxt.log.write("\t-> File not compressed... Skipping decompression of file {}".format(f.fullpath))
				
		for p in decompress_procs:
			p.communicate()
	#end run()
#end class GeneratePseudoreplicates