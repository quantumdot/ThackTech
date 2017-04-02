import os
from ThackTech.Pipelines import PipelineModule


class RemoveDecompressedFiles(PipelineModule):
	
	def __init__(self, **kwargs):
		super(RemoveDecompressedFiles, self).__init__('RemoveDecompressedFiles', 'Remove Decompressed Files', **kwargs)
	#end __init__()
	
	def run(self, cxt):
		for file_group in cxt.sample.get_attribute('decompressed_files'):
			decompressed_file = cxt.sample.find_files(lambda f: f.fullpath == file_group[0])[0]
			decompressed_file._set_path(file_group[1])
			cxt.log.write("\t-> Removing decompressed file {}....\n".format(file_group[0]))
			cxt.log.flush()
			os.remove(file_group[0])
	#end run()
#end class GeneratePseudoreplicates