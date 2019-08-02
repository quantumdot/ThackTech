import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class BamCoverage(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='BamCoverage', short_description='Convert Bam to coverage track')
		super_args.update(**kwargs)
		super(BamCoverage, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('bamCoverage_path', str, 'bamCoverage', desc="Path to the bamCoverage program"))
		
		#output options
		self.add_parameter(ModuleParameter('output_format', str, 'bigwig', choices=['bigwig', 'bedgraph'], desc="Output file format."))
		
		#Optional arguments
		self.add_parameter(ModuleParameter('scale_factor', float, 1.0, desc="Scalar to muliply against the computed scaling factor"))
		self.add_parameter(ModuleParameter('MNase', bool, False, desc="process as MNase data"))
		self.add_parameter(ModuleParameter('offset', int, None, nullable=True))
		self.add_parameter(ModuleParameter('filter_strand', str, None, nullable=True, choices=['forward', 'reverse'], desc=""))
		self.add_parameter(ModuleParameter('region', str, None, nullable=True))
		self.add_parameter(ModuleParameter('blacklist', str, None, nullable=True))
		self.add_parameter(ModuleParameter('bin_size', int, 1, desc="Size of the bins, in bases, for the output."))
			
		
		#Read coverage normalization options:
		self.add_parameter(ModuleParameter('normalize', str, 'RPKM', nullable=True, choices=['RPKM','CPM','BPM','RPGC'], desc="Normalization type to perform"))
		self.add_parameter(ModuleParameter('ignore_contigs', list, None, nullable=True, desc="list of contig (chromosome) names to exclude for normalization computation."))
		self.add_parameter(ModuleParameter('skip_na', bool, False, desc="This parameter determines if non-covered regions (regions without overlapping reads) in a BAM file should be skipped."))
		self.add_parameter(ModuleParameter('smooth', int, None, nullable=True, desc="The smooth length defines a window, larger than the binSize, to average the number of reads."))
		
		#Read processing options:
		self.add_parameter(ModuleParameter('extend_reads', int, None, nullable=True))
		self.add_parameter(ModuleParameter('ignore_duplicates', bool, False))
		self.add_parameter(ModuleParameter('min_quality', int, None, nullable=True))
		self.add_parameter(ModuleParameter('center_reads', bool, False))
		self.add_parameter(ModuleParameter('sam_flag_include', int, None, nullable=True))
		self.add_parameter(ModuleParameter('sam_flag_exclude', int, None, nullable=True))
		self.add_parameter(ModuleParameter('min_frag_len', int, None, nullable=True))
		self.add_parameter(ModuleParameter('max_frag_len', int, None, nullable=True))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bam')
	#end __declare_resolvers()

	def tool_versions(self):
		bcp = self.get_parameter_value_as_string('bamCoverage_path')
		return {
			'bamCoverage': self._call_output(bcp+r" --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		bam = self.resolve_input('bam', cxt)
		if bam is None or (isinstance(bam, list) and len(bam) <= 0):
			cxt.log.write('No file returned for resolver "bam", exiting...\n')
			cxt.log.flush()
			return
		
		if isinstance(bam, list):
			cxt.log.write('WARNING: resolver returned a list with {} item(s), but this module only supports one input from this resolver. Taking the first item.\n'.format(len(bam)))
			bam = bam[0]
		
		out_ext = ('bg' if self.get_parameter_value_as_string('output_format') == 'bedgraph' else 'bw')
		sample_basename = bam.basename_with_ext('rpkm_norm.{}'.format(out_ext))
		bamcoverage_args = [
			self.get_parameter_value_as_string('bamCoverage_path'),
			'--bam', bam.fullpath,
			'--outFileName', os.path.join(cxt.sample.dest, sample_basename),
			'--outFileFormat', self.get_parameter_value_as_string('output_format'),
			'--binSize', self.get_parameter_value_as_string('bin_size'),
			'--numberOfProcessors', str(self.processors),
			'--scaleFactor', self.get_parameter_value_as_string('scale_factor'),
			#'--verbose'
		]
		if self.get_parameter_value('normalize') is not None:
			bamcoverage_args.extend(['--normalizeUsing', self.get_parameter_value_as_string('normalize')])

		if self.get_parameter_value('MNase'):
			bamcoverage_args.append('--MNase')

		if self.get_parameter_value('filter_strand') is not None:
			bamcoverage_args.extend(['--filterRNAstrand', self.get_parameter_value_as_string('filter_strand')])

		if self.get_parameter_value('smooth') is not None:
			bamcoverage_args.extend(['--smoothLength', self.get_parameter_value_as_string('smooth')])

		if self.get_parameter_value('skip_na'):
			bamcoverage_args.append('--skipNonCoveredRegions')

		cxt.log.write('Generating normalized signal track....')
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bamcoverage_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bamcoverage_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return {
			'coverage_signal': os.path.join(cxt.sample.dest, sample_basename)
		}
	#end run()
#end class RPKMNormBigWig