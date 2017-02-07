import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class RemoveDecompressedFiles(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'RemoveDecompressedFiles', 'Remove Decompressed Files')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		for f in sample.get_attribute('decompressed_files'):
			logfile.write("\t-> Removing decompressed file %s....\n" % (f,))
			logfile.flush()
			os.remove(f)
	#end run()
#end class GeneratePseudoreplicates