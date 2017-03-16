import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class GenerateMACsModelFigure(PipelineModule):
	def __init__(self):
		PipelineModule.__init__(self, 'MacsModelFig', 'Generate MACS Model Figure')
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):
		model_path = os.path.join(sample.dest, sample.name+'_model.r')
		if os.path.exists(model_path):
			logfile.write('Generating MACS model figure.....\n')
			logfile.flush()
			with open(model_path, 'r') as rmodel:
				with open(os.devnull, 'w') as devnull:
					self._run_subprocess(['R', '--quiet', '--vanilla'], stdin=rmodel, stderr=subprocess.STDOUT, stdout=devnull, cwd=sample.dest)
			return { 'model_figure': os.path.join(sample.dest, sample.name+'_model.pdf') }
		else:
			logfile.write('No model figure to generate.....\n')
			logfile.flush()
			return None
	#end run()
#end class GenerateMACsModelFigure