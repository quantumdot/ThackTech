import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class BamToRpkmNormBigWig(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'BamToRpkmNormBigWig', 'Make RPKM Norm BigWig from BAM')
		
		self.add_parameter(ModuleParameter('output_format', str, 	'bigwig',	desc="Output file type. Either 'bigwig' or 'bedgraph'."))
		self.add_parameter(ModuleParameter('bin_size', 		int, 	1,			desc="Size of the bins, in bases, for the output of the bigwig/bedgraph file."))
		
		self._name_resolver('bam')
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	def run(self, sample, logfile):
		#dest = sample.get_attribute('origional_dest') if sample.has_attribute('origional_dest') else sample.dest
		sample_basename = "%s.%s.%s" % (sample.name, 'rpkm_norm', ('bg' if self.get_parameter_value_as_string('output_format') == 'bedgraph' else 'bw'))
		bamcoverage_args = [
			'bamCoverage',
			'--bam', self.resolve_input('bam', sample),
			'--outFileName', os.path.join(sample.dest, sample_basename),
			'--outFileFormat', self.get_parameter_value_as_string('output_format'),
			'--binSize', self.get_parameter_value_as_string('bin_size'),
			'--normalizeUsingRPKM',
			'--numberOfProcessors', str(self.processors),
			'--verbose'
		]
		logfile.write('Generating normalized signal track....')
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(bamcoverage_args))
		logfile.write("\n..............................................\n")
		logfile.flush()
		self._run_subprocess(bamcoverage_args, cwd=sample.dest, stderr=subprocess.STDOUT, stdout=logfile)
		
		return {
			'rpkm_norm_bw': os.path.join(sample.dest, sample_basename)
		}
	#end run()
#end class RPKMNormBigWig