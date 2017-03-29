from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class DecompressFiles(PipelineModule):
	
	def __init__(self, **kwargs):
		super(DecompressFiles, self).__init__('DecompressFiles', 'Decompress Files', **kwargs)
		self._name_resolver('files')
		self.add_parameter(ModuleParameter('overwrite',	bool, False, desc="Overwrite the decompressed file if it already exists."))
	#end __init__()

	
	def run(self, cxt):
		cxt.sample.add_attribute('decompressed_files', [])
		file_lookups = self.resolve_input('files', cxt)
		files = []
		decompress_procs = []
		for fl in file_lookups:
			ofiles = cxt.sample.get_file(fl[0], fl[1])
			if isinstance(ofiles, str):
				#just a single file here
				if filetools.is_compressed(ofiles):
					cxt.log.write("\t-> Decompressing file %s....\n" % (ofiles,))
					cxt.log.flush()
					d, p = filetools.extract(ofiles, cxt.sample.dest, overwrite=self.get_parameter_value('overwrite'))
					cxt.sample.add_file(fl[0], fl[1], d)
					cxt.sample.add_attribute('decompressed_files', cxt.sample.get_attribute('decompressed_files') + [d])
					decompress_procs.append(p)
			else:
				dfiles = []
				for of in ofiles:
					if filetools.is_compressed(of):
						cxt.log.write("\t-> Decompressing file %s....\n" % (of,))
						cxt.log.flush()
						d, p = filetools.extract(of, cxt.sample.dest, overwrite=self.get_parameter_value('overwrite'))
						dfiles.append(d)
						cxt.sample.add_attribute('decompressed_files', cxt.sample.get_attribute('decompressed_files') + [d])
						decompress_procs.append(p)
				cxt.sample.add_file(fl[0], fl[1], dfiles)
		for p in decompress_procs:
			p.communicate()
		return files
	#end run()
#end class GeneratePseudoreplicates