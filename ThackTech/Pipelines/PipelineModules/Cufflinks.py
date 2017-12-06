import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class Cufflinks(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Cufflinks', short_description='Transcript assembly with Cufflinks')
		super_args.update(**kwargs)
		super(Cufflinks, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		cl_path = self.get_parameter_value_as_string('cufflinks_path')
		return {
			'cufflinks': self._call_output(cl_path+" 2>&1 | perl -ne 'if(m/.*cufflinks.*v([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('cufflinks_path', str, 'cufflinks', desc="Path to the cufflinks executable."))
		
		#general options
		ltype_opts = ['fr-unstranded', 'ff-firststrand', 'ff-secondstrand', 'ff-unstranded', 'fr-firststrand', 'fr-secondstrand', 'transfrags']
		self.add_parameter(ModuleParameter('library_type', str, ltype_opts[0], choices=ltype_opts, desc="Sets the library type. See --library-type from cufflinks."))
		lnorm_opts = ['classic-fpkm']
		self.add_parameter(ModuleParameter('library_norm', str, lnorm_opts[0], choices=lnorm_opts, desc="Method used to normalize library sizes. See --library-norm-method from cufflinks."))
		
		# @todo: the following are not currently implemented!
		#   --seed                       value of random number generator seed                 [ default:      0 ]
		#   -G/--GTF                     quantitate against reference transcript annotations
		#   -g/--GTF-guide               use reference transcript annotation to guide assembly
		#   -M/--mask-file               ignore all alignment within transcripts in this file
		#   -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
		#   -u/--multi-read-correct      use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]

		
		#Advanced Abundance Estimation Options:
		self.add_parameter(ModuleParameter('mean_frag_len', int, 200, desc="average fragment length (unpaired reads only). See -m from cufflinks."))
		self.add_parameter(ModuleParameter('std_frag_len', int, 80, desc="fragment length std deviation (unpaired reads only). See -s from cufflinks."))
		self.add_parameter(ModuleParameter('max_mle_iterations', int, 5000, desc="maximum iterations allowed for MLE calculation. See --max-mle-iterations from cufflinks"))
		self.add_parameter(ModuleParameter('hits_norm', str, 'total', choices=['total', 'compatible'], desc="Hit normalization method. See --total-hits-norm and --compatible-hits-norm from cufflinks."))
		self.add_parameter(ModuleParameter('num_frag_count_draws', int, 100, desc="Number of fragment generation samples"))
		self.add_parameter(ModuleParameter('num_frag_assign_draws', int, 50, desc="Number of fragment assignment samples per generation"))
		self.add_parameter(ModuleParameter('max_frag_multihits', int, None, nullable=True, desc="Maximum number of alignments allowed per fragment, None for unlimited"))
		self.add_parameter(ModuleParameter('no_effective_length_correction', bool, False, desc="No effective length correction"))
		self.add_parameter(ModuleParameter('no_length_correction', bool, False, desc="No length correction"))

		#Advanced Assembly Options:
		self.add_parameter(ModuleParameter('label', str, 'CUFF', desc="assembled transcripts have this ID prefix"))
		self.add_parameter(ModuleParameter('min_isoform_fraction', float, 0.10, desc="suppress transcripts below this abundance level"))
		self.add_parameter(ModuleParameter('pre_mrna_fraction', float, 0.15, desc="suppress intra_intronic transcripts below this level"))
		self.add_parameter(ModuleParameter('max_intron_length', int, 300000, desc="ignore alignments with gaps longer than this"))
		self.add_parameter(ModuleParameter('junc_alpha', float, 0.001, desc="alpha for junction binomial test filter"))
		self.add_parameter(ModuleParameter('small_anchor_fraction', float, 0.09, desc="percent read overhang taken as 'suspiciously small'"))
		self.add_parameter(ModuleParameter('min_frags_per_transfrag', int, 10, desc="minimum number of fragments needed for new transfrags"))
		self.add_parameter(ModuleParameter('overhang_tolerance', int, 8, desc="number of terminal exon bp to tolerate in introns"))
		self.add_parameter(ModuleParameter('max_bundle_length', int, 3500000, desc="maximum genomic length allowed for a given bundle"))
		self.add_parameter(ModuleParameter('max_bundle_frags', int, 500000, desc="maximum fragments allowed in a bundle before skipping"))
		self.add_parameter(ModuleParameter('min_intron_length', int, 50, desc="minimum intron size allowed in genome"))
		self.add_parameter(ModuleParameter('trim_3_avgcov_thresh', int, 10, desc="minimum avg coverage required to attempt 3' trimming"))
		self.add_parameter(ModuleParameter('trim_3_dropoff_frac', float, 0.1, desc="fraction of avg coverage below which to trim 3' end"))
		self.add_parameter(ModuleParameter('max_multiread_fraction', float, 0.75, desc="maximum fraction of allowed multireads per transcript"))
		self.add_parameter(ModuleParameter('overlap_radius', int, 50, desc="maximum gap size to fill between transfrags (in bp)"))
		
		#Advanced Reference Annotation Guided Assembly Options:
		self.add_parameter(ModuleParameter('no_faux_reads', bool, False, desc="disable tiling by faux reads"))
		self.add_parameter(ModuleParameter('3_overhang_tolerance', int, 600, desc="overhang allowed on 3' end when merging with reference"))
		self.add_parameter(ModuleParameter('intron_overhang_tolerance', int, 30, desc="overhang allowed inside reference intron when merging"))
		
		
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
		self._name_resolver('guide_gff')
	#end __declare_resolvers()
	
	def run(self, cxt):
		cl_args = [
			self.get_parameter_value_as_string('cufflinks_path'),
			'--quiet',
			'--num-threads', str(self.processors),
			'--output-dir', cxt.sample.dest,
			
			'--library-type', self.get_parameter_value_as_string('library_type'),
			'--library-norm-method', self.get_parameter_value_as_string('library_norm'),
			
			'--max-mle-iterations', self.get_parameter_value_as_string('max_mle_iterations'),
			'--{}-hits-norm'.format(self.get_parameter_value_as_string('hits_norm')),
			'--num-frag-count-draws', self.get_parameter_value_as_string('num_frag_count_draws'),
			'--num-frag-assign-draws', self.get_parameter_value_as_string('num_frag_assign_draws'),
			
			'--label', self.get_parameter_value_as_string('label'),
			'--min-isoform-fraction', self.get_parameter_value_as_string('min_isoform_fraction'),
			'--pre-mrna-fraction', self.get_parameter_value_as_string('pre_mrna_fraction'),
			'--max-intron-length', self.get_parameter_value_as_string('max_intron_length'),
			'--junc-alpha', self.get_parameter_value_as_string('junc_alpha'),
			'--small-anchor-fraction', self.get_parameter_value_as_string('small_anchor_fraction'),
			'--min-frags-per-transfrag', self.get_parameter_value_as_string('min_frags_per_transfrag'),
			'--overhang-tolerance', self.get_parameter_value_as_string('overhang_tolerance'),
			'--max-bundle-length', self.get_parameter_value_as_string('max_bundle_length'),
			'--max-bundle-frags', self.get_parameter_value_as_string('max_bundle_frags'),
			'--min-intron-length', self.get_parameter_value_as_string('min_intron_length'),
			'--trim-3-avgcov-thresh', self.get_parameter_value_as_string('trim_3_avgcov_thresh'),
			'--trim-3-dropoff-frac', self.get_parameter_value_as_string('trim_3_dropoff_frac'),
			'--max-multiread-fraction', self.get_parameter_value_as_string('max_multiread_fraction'),
			'--overlap-radius', self.get_parameter_value_as_string('overlap_radius'),
			
			'--3-overhang-tolerance', self.get_parameter_value_as_string('3_overhang_tolerance'),
			'--intron-overhang-tolerance', self.get_parameter_value_as_string('intron_overhang_tolerance'),
		]
		
		if self.get_parameter_value('no_effective_length_correction'):
			cl_args.append('--no-effective-length-correction')
			
		if self.get_parameter_value('no_length_correction'):
			cl_args.append('--no-length-correction')
			
		if self.get_parameter_value('no_faux_reads'):
			cl_args.append('--no-faux-reads')
		
		if self.get_parameter_value('max_frag_multihits') is not None:
			cl_args.extend(['--max-frag-multihits', self.get_parameter_value_as_string('max_frag_multihits')])
		
		if not cxt.sample.get_attribute('PE'):
			cl_args.extend([
				'–-frag-len-mean', self.get_parameter_value_as_string('mean_frag_len'),
				'--–frag-len-std-dev', self.get_parameter_value_as_string('std_frag_len')
			])
			
		#add the input alignments that we are quantifying/assembling
		alignments = self.resolve_input('alignments', cxt)
		for a in alignments:
			cl_args.append(a.fullpath)
			break #only support one!
		
		
		#OK, we now have all the arguments setup, lets actually run cufflinks
		cxt.log.write("\t-> Assembling transcripts with cufflinks......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(cl_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(cl_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		
		return [
			FileInfo(os.path.join(cxt.sample.dest, 'transcripts.gtf'), FileContext.from_module_context(cxt, 'assembled_transcripts')),
			FileInfo(os.path.join(cxt.sample.dest, 'isoforms.fpkm_tracking'), FileContext.from_module_context(cxt, 'transcript_fpkm_tracking')),
			FileInfo(os.path.join(cxt.sample.dest, 'genes.fpkm_tracking'), FileContext.from_module_context(cxt, 'gene_fpkm_tracking'))
		]
		
	#end run()
#end class HelloWorld