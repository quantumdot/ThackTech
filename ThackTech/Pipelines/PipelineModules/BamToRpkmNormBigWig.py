import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class BamToRpkmNormBigWig(PipelineModule):
	
	def __init__(self, **kwargs):
		super(BamToRpkmNormBigWig, self).__init__('BamToRpkmNormBigWig', 'Make RPKM Norm BigWig from BAM', **kwargs)
		
		self.add_parameter(ModuleParameter('output_format', str, 	'bigwig',	desc="Output file type. Either 'bigwig' or 'bedgraph'."))
		self.add_parameter(ModuleParameter('bin_size', 		int, 	1,			desc="Size of the bins, in bases, for the output of the bigwig/bedgraph file."))
		
		self._name_resolver('bam')
	#end __init__()

	def tool_versions(self):
		return {
			'bamCoverage': self._call_output("bamCoverage --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		bam = self.resolve_input('bam', cxt)
		out_ext = ('bg' if self.get_parameter_value_as_string('output_format') == 'bedgraph' else 'bw')
		#dest = cxt.sample.get_attribute('origional_dest') if cxt.sample.has_attribute('origional_dest') else cxt.sample.dest
		sample_basename = bam.basename_with_ext('rpkm_norm.{}'.format(out_ext))
		bamcoverage_args = [
			'bamCoverage',
			'--bam', bam.fullpath,
			'--outFileName', os.path.join(cxt.sample.dest, sample_basename),
			'--outFileFormat', self.get_parameter_value_as_string('output_format'),
			'--binSize', self.get_parameter_value_as_string('bin_size'),
			'--normalizeUsingRPKM',
			'--numberOfProcessors', str(self.processors),
			#'--verbose'
		]
		cxt.log.write('Generating normalized signal track....')
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bamcoverage_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bamcoverage_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return {
			'rpkm_norm_bw': os.path.join(cxt.sample.dest, sample_basename)
		}
	#end run()
#end class RPKMNormBigWig