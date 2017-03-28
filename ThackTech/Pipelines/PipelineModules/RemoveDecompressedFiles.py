import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class RemoveDecompressedFiles(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'RemoveDecompressedFiles', 'Remove Decompressed Files')
	#end __init__()
	
	def run(self, cxt):
		for f in cxt.sample.get_attribute('decompressed_files'):
			cxt.log.write("\t-> Removing decompressed file %s....\n" % (f,))
			cxt.log.flush()
			os.remove(f)
	#end run()
#end class GeneratePseudoreplicates