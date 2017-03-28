import os
import random
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter
import pysam


class GeneratePseudoreplicates(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'Pseudoreplicates', 'Generate Pseudoreplicates')
		self._name_resolver('alignments')
		self.add_parameter(ModuleParameter('nreplicates', int, 2))
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
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