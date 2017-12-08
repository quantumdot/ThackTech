import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class Cuffquant(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cuffquant', short_description='Quantify transcript expression with Cuffquant')
		super_args.update(**kwargs)
		super(Cuffquant, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		cq_path = self.get_parameter_value_as_string('cuffquant_path')
		return {
			'cuffquant': self._call_output(cq_path+" 2>&1 | perl -ne 'if(m/.*cuffquant.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('cuffquant_path', str, 'cuffquant', desc="Path to the cuffquant executable."))
		
		#General Options:
		ltype_opts = ['fr-unstranded', 'ff-firststrand', 'ff-secondstrand', 'ff-unstranded', 'fr-firststrand', 'fr-secondstrand', 'transfrags']
		self.add_parameter(ModuleParameter('library_type', str, ltype_opts[0], choices=ltype_opts, desc="Sets the library type. See --library-type from cuffquant."))
		self.add_parameter(ModuleParameter('use_mask_file', bool, False, desc="ignore all alignment within transcripts in 'mask' file. Requires that the 'mask' resolver is set."))
		self.add_parameter(ModuleParameter('frag_bias_correct', bool, False, desc="use bias correction - reference fasta required."))
		self.add_parameter(ModuleParameter('multi_read_correct', bool, False, desc="use 'rescue method' for multi-reads"))
		
		#Advanced Options:
		self.add_parameter(ModuleParameter('mean_frag_len', int, 200, desc="average fragment length (unpaired reads only). See -m from cuffquant."))
		self.add_parameter(ModuleParameter('std_frag_len', int, 80, desc="fragment length std deviation (unpaired reads only). See -s from cuffquant."))
		self.add_parameter(ModuleParameter('min_alignment_count', int, 10, desc="minimum number of alignments in a locus for testing."))
		self.add_parameter(ModuleParameter('max_mle_iterations', int, 5000, desc="maximum iterations allowed for MLE calculation. See --max-mle-iterations from cuffquant"))
		self.add_parameter(ModuleParameter('max_bundle_frags', int, 500000, desc="maximum fragments allowed in a bundle before skipping"))
		self.add_parameter(ModuleParameter('max_frag_multihits', int, None, nullable=True, desc="Maximum number of alignments allowed per fragment, None for unlimited"))
		self.add_parameter(ModuleParameter('no_effective_length_correction', bool, False, desc="No effective length correction"))
		self.add_parameter(ModuleParameter('no_length_correction', bool, False, desc="No length correction"))

	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('mask')
		self._name_resolver('ref_fasta')
		self._name_resolver('alignments')
		self._name_resolver('merged_gff')
	#end __declare_resolvers()
	
	def run(self, cxt):
		cq_args = [
			self.get_parameter_value_as_string('cuffquant_path'),
			'--quiet', '--no-update-check',
			'--num-threads', str(self.processors),
			'--output-dir', cxt.sample.dest,
			
			'--min-alignment-count', self.get_parameter_value_as_string('min_alignment_count'),
			'--max-mle-iterations', self.get_parameter_value_as_string('max_mle_iterations'),
			'--max-bundle-frags', self.get_parameter_value_as_string('max_bundle_frags'),
		]
		
		if self.get_parameter_value('use_mask_file'):
			mf = self.resolve_input('mask', cxt)
			cq_args.extend(['--mask-file', mf])
			
		if self.get_parameter_value('frag_bias_correct'):
			rf = self.resolve_input('ref_fasta', cxt)
			cq_args.extend(['--frag-bias-correct', rf])
			
		if self.get_parameter_value('multi_read_correct'):
			cq_args.append('--multi-read-correct')
			
		if self.get_parameter_value('max_frag_multihits') is not None:
			cq_args.extend(['--max-frag-multihits', self.get_parameter_value_as_string('max_frag_multihits')])
		
		if not cxt.sample.get_attribute('PE'):
			cq_args.extend([
				'–-frag-len-mean', self.get_parameter_value_as_string('mean_frag_len'),
				'--–frag-len-std-dev', self.get_parameter_value_as_string('std_frag_len')
			])
		
		if self.get_parameter_value('no_effective_length_correction'):
			cq_args.append('--no-effective-length-correction')
			
		if self.get_parameter_value('no_length_correction'):
			cq_args.append('--no-length-correction')
		
		#add the merged gtf that we are quantifying/assembling
		gtfs = self.resolve_input('merged_gff', cxt)
		for g in gtfs:
			cq_args.append(g.fullpath)
			break #only support one!
			
		#add the input alignments that we are quantifying/assembling
		alignments = self.resolve_input('alignments', cxt)
		for a in alignments:
			cq_args.append(a.fullpath)
			break #only support one!
		
		#OK, we now have all the arguments setup, lets actually run cuffquant
		cxt.log.write("\t-> Quantifying gene and transcript expression with cuffquant......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(cq_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(cq_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		
		return [
			FileInfo(os.path.join(cxt.sample.dest, 'abundances.cxb'), FileContext.from_module_context(cxt, 'quantified_transcripts'))
		]
	#end run()
#end class HelloWorld