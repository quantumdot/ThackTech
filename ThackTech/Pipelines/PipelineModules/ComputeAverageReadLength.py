import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ComputeAverageReadLength(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='ComputeAverageReadLength', short_description='Compute Average Read Length')
		super_args.update(**kwargs)
		super(ComputeAverageReadLength, self).__init__(**super_args)
		
		self._name_resolver('fastq')
		
		self.add_parameter(ModuleParameter('script_path', str, '/home/josh/scripts/avg_read_len_fq.awk', desc="Location of average read length script."))
	#end __init__()

	
	def run(self, cxt):	
		read_length = subprocess.check_output(['awk', '-f', self.get_parameter_value_as_string('script_path'), self.resolve_input('fastq', cxt)])
		cxt.log.write("-> Average Read Length = %d\n" % (int(read_length),))
		cxt.log.flush()
		cxt.cxt.sample.set_attribute('avgreadlength', int(read_length))
	#end run()
#end class GeneratePseudoreplicates