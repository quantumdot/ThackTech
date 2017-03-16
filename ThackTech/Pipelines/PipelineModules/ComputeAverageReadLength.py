import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ComputeAverageReadLength(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'ComputeAverageReadLength', 'Compute Average Read Length')
		self._name_resolver('fastq')
		self.add_parameter(ModuleParameter('script_path', str, '/home/josh/scripts/avg_read_len_fq.awk', desc="Location of average read length script."))
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, sample, logfile):	
		read_length = subprocess.check_output(['awk', '-f', self.get_parameter_value_as_string('script_path'), self.resolve_input('fastq', sample)])
		logfile.write("-> Average Read Length = %d\n" % (int(read_length),))
		logfile.flush()
		sample.add_attribute('avgreadlength', int(read_length))
	#end run()
#end class GeneratePseudoreplicates