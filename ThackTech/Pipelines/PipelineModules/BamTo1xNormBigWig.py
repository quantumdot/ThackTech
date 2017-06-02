import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class RPKMNormBigWig(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='RPKMBigwig', short_description='RPKM Norm BigWig')
		super_args.update(**kwargs)
		super(RPKMNormBigWig, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('output_format', str, 	'bigwig',	desc="Output file type. Either 'bigwig' or 'bedgraph'."))
		self.add_parameter(ModuleParameter('bin_size', 		int, 	1,			desc="Size of the bins, in bases, for the output of the bigwig/bedgraph file."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bam')
	#end __declare_resolvers()

	def tool_versions(self):
		return {
			'bamCoverage': self._call_output("bamCoverage --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		#dest = cxt.sample.get_attribute('origional_dest') if cxt.sample.has_attribute('origional_dest') else cxt.sample.dest
		cxt.sample_basename = "%s.%s.%s" % (cxt.sample.name, 'rpkm_norm', ('bg' if self.get_parameter_value_as_string('output_format') == 'bedgraph' else 'bw'))
		bamcoverage_args = [
			'bamCoverage',
			'--bam', self.resolve_input('bam', cxt),
			'--outFileName', os.path.join(cxt.sample.dest, cxt.sample_basename),
			'--outFileFormat', self.get_parameter_value_as_string('output_format'),
			'--binSize', self.get_parameter_value_as_string('bin_size'),
			'--normalizeUsingRPKM',
			'--numberOfProcessors', str(self.processors),
			'--verbose'
		]
		cxt.log.write('Generating normalized signal track....')
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bamcoverage_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bamcoverage_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return {
			'rpkm_norm_bw': os.path.join(cxt.sample.dest, cxt.sample_basename)
		}
	#end run()
#end class RPKMNormBigWig