import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class DecompressFiles(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'DecompressFiles', 'Decompress Files')
		self._name_resolver('files')
		self.add_parameter(ModuleParameter('overwrite',	bool, False, desc="Overwrite the decompressed file if it already exists."))
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		sample.add_attribute('decompressed_files', [])
		file_lookups = self.resolve_input('files', sample)
		files = []
		decompress_procs = []
		for fl in file_lookups:
			ofiles = sample.get_file(fl[0], fl[1])
			if isinstance(ofiles, str):
				#just a single file here
				if filetools.is_compressed(ofiles):
					logfile.write("\t-> Decompressing file %s....\n" % (ofiles,))
					logfile.flush()
					d, p = filetools.extract(ofiles, sample.dest, overwrite=self.get_parameter_value('overwrite'))
					sample.add_file(fl[0], fl[1], d)
					sample.add_attribute('decompressed_files', sample.get_attribute('decompressed_files') + [d])
					decompress_procs.append(p)
			else:
				dfiles = []
				for of in ofiles:
					if filetools.is_compressed(of):
						logfile.write("\t-> Decompressing file %s....\n" % (of,))
						logfile.flush()
						d, p = filetools.extract(of, sample.dest, overwrite=self.get_parameter_value('overwrite'))
						dfiles.append(d)
						sample.add_attribute('decompressed_files', sample.get_attribute('decompressed_files') + [d])
						decompress_procs.append(p)
				sample.add_file(fl[0], fl[1], dfiles)
		for p in decompress_procs:
			p.communicate()
		return files
	#end run()
#end class GeneratePseudoreplicates