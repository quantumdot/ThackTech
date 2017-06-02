import os
import random
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter
import pysam


class GeneratePseudoreplicates(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Pseudoreplicates', short_description='Generate Pseudoreplicates')
		super_args.update(**kwargs)
		super(GeneratePseudoreplicates, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('nreplicates', int, 2))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
	#end __declare_resolvers()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'pseudoreps')
		filetools.ensure_dir(dest_dir)
		alignments = self.resolve_input('alignments', cxt.sample)
		samfile = pysam.AlignmentFile(alignments, 'rb')
		
		replicates = []
		rep_count = self.get_parameter_value('nreplicates')
		out_files = {}
		for i in range(rep_count):
			key = 'pr' + str(i+1)
			out_files[key] = os.path.join(dest_dir, cxt.sample.name+'.'+key+'.bam')
			replicates.append(pysam.AlignmentFile(out_files[key], 'wb', template=samfile))
		
		for read in samfile.fetch():
			replicates[random.randint(0, rep_count-1)].write(read)
		
		samfile.close()
		for rep in replicates:
			rep.close()

		return out_files
	#end run()
#end class GeneratePseudoreplicates