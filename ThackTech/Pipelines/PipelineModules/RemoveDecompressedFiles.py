import os
from ThackTech.Pipelines import PipelineModule


class RemoveDecompressedFiles(PipelineModule):
	
	def __init__(self, **kwargs):
		super(RemoveDecompressedFiles, self).__init__('RemoveDecompressedFiles', 'Remove Decompressed Files', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		for f in cxt.sample.get_attribute('decompressed_files'):
			cxt.log.write("\t-> Removing decompressed file %s....\n" % (f[],))
			cxt.log.flush()
			os.remove(f)
	#end run()
#end class GeneratePseudoreplicates