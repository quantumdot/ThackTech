from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class MoveFiles(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MoveFiles', short_description='Move Files')
		super_args.update(**kwargs)
		super(MoveFiles, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('files')
		self._name_resolver('dest_dir')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {}
	#end tool_versions()
	
	def run(self, cxt):
		outdir = self.resolve_input_path('dest_dir', cxt)
		filetools.ensure_dir(outdir)
		
		files_to_move = self.resolve_input('files', cxt)
		for f in files_to_move:
			cxt.log.write('Moving file: {} -> {}'.format(f.fullpath, outdir))
			f.move(outdir)
	#end run()
#end class MoveFiles	