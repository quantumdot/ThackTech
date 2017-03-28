import os
import subprocess
from ThackTech.Pipelines import PipelineModule


class GenerateMACsModelFigure(PipelineModule):
	def __init__(self, **kwargs):
		super(GenerateMACsModelFigure, self).__init__('MacsModelFig', 'Generate MACS Model Figure', **kwargs)
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		model_path = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')
		if os.path.exists(model_path):
			cxt.log.write('Generating MACS model figure.....\n')
			cxt.log.flush()
			with open(model_path, 'r') as rmodel:
				with open(os.devnull, 'w') as devnull:
					self._run_subprocess(['R', '--quiet', '--vanilla'], stdin=rmodel, stderr=subprocess.STDOUT, stdout=devnull, cwd=cxt.sample.dest)
			return { 'model_figure': os.path.join(cxt.sample.dest, cxt.sample.name+'_model.pdf') }
		else:
			cxt.log.write('No model figure to generate.....\n')
			cxt.log.flush()
			return None
	#end run()
#end class GenerateMACsModelFigure