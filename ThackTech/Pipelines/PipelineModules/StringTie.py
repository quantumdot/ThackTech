import os
import platform
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class StringTieBase(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='StringTie', short_description='Transcript assembly with StringTie.')
		super_args.update(**kwargs)
		super(StringTieBase, self).__init__(**super_args)
	#end __init__()
	
	def tool_versions(self):
		return {
			'stringtie': self._call_output([self.get_parameter_value_as_string('stringtie_path'), '--version'], stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('stringtie_path', str, 'stringtie', desc="Path to the StringTie executable."))
	#end __declare_parameters()

#end class StringTieBase


class StringTieQuant(StringTieBase):
	def __init__(self, **kwargs):
		super_args = dict(name='StringTieQuant', short_description='Assemble transcripts with StringTie.')
		super_args.update(**kwargs)
		super(StringTieQuant, self).__init__(**super_args)
	#end __init__()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
		self._name_resolver('guide_gff')
	#end __declare_resolvers()

	def _declare_parameters(self):
		super(StringTieMerge, self)._declare_parameters()
		
		#--fr/--rf
		self.add_parameter(ModuleParameter('stranded', str, None, nullable=True, choices=['rf', 'fr'], desc="Sets the library type. See --rf and --fr from stringtie."))
		
		#FLAGS
		#-v
		self.add_parameter(ModuleParameter('verbose', bool, False, desc="log bundle processing details. See -v from stringtie"))
		#-t
		self.add_parameter(ModuleParameter('cov_trim', bool, True, desc="trimming of predicted transcripts based on coverage. See -t from stringtie"))
		#-u
		self.add_parameter(ModuleParameter('mm_corr', bool, True, desc="multi-mapping correction. See -u from stringtie"))
		#-e
		self.add_parameter(ModuleParameter('estimate', bool, False, desc="Limits processing for read alignments to only estimate quantity. See -e from stringtie."))
		
		
		#Value-based options
		#-f
		self.add_parameter(ModuleParameter('min_abd', float, 0.1, desc="minimum isoform fraction. See -f from stringtie."))
		#-m
		self.add_parameter(ModuleParameter('min_len', int, 200, desc="minimum assembled transcript length. See -m from stringtie."))
		#-a
		self.add_parameter(ModuleParameter('min_anc', int, 10, desc="minimum anchor length for junctions. See -a from stringtie."))
		#-j
		self.add_parameter(ModuleParameter('min_jun', int, 1, desc="minimum junction coverage. See -j from stringtie."))
		#-c
		self.add_parameter(ModuleParameter('min_cov', float, 2.5, desc="minimum reads per bp coverage to consider for transcript assembly. See -c from stringtie."))
		#-g
		self.add_parameter(ModuleParameter('bun_gap', int, 50, desc="gap between read mappings triggering a new bundle. See -g from stringtie."))
		#-M
		self.add_parameter(ModuleParameter('max_mhf', float, 0.95, desc="fraction of bundle allowed to be covered by multi-hit reads. See -M from stringtie."))
		#-x
		self.add_parameter(ModuleParameter('skip_scaffold', list, None, nullable=True, desc="List of scaffold names (e.x. chrM) to ignore for processing. See -x from stringtie."))
		
		
		#output types
		#-B
		self.add_parameter(ModuleParameter('out_ballgown', bool, False, desc="Output Ballgown compatible files. See -B from stringtie"))
		#-A
		self.add_parameter(ModuleParameter('out_abundance', bool, False, desc="Output gene abundance data. See -A from stringtie"))
		#-C
		self.add_parameter(ModuleParameter('out_full_cov', bool, False, desc="Output reference transcripts that are fully covered. See -C from stringtie"))
		
	#end _declare_parameters()
	
	
	def run(self, cxt):
		gtf_out = os.path.join(cxt.sample.dest, cxt.sample.name+'.gtf')
		output_files = [
			FileInfo(gtf_out, FileContext.from_module_context(cxt, 'assembled_transcripts'))
		]	
		
		
		st_args = [
			self.get_parameter_value_as_string('stringtie_path'),
			'-p', str(self.processors), #not sure if this is actually supported, but lets give it a shot!
			'-l', cxt.sample.name,
			'-o', gtf_out,
			#'-G', #<guide_gff>, # @todo: need to implement this!
			
			#optional options
			'-f', self.get_parameter_value_as_string('min_abd'),
			'-m', self.get_parameter_value_as_string('min_len'),
			'-a', self.get_parameter_value_as_string('min_anc'),
			'-j', self.get_parameter_value_as_string('min_jun'),
			'-c', self.get_parameter_value_as_string('min_cov'),
			'-g', self.get_parameter_value_as_string('bun_gap'),
			'-M', self.get_parameter_value_as_string('max_mhf'),
		]
		guide_gff = self.resolve_input('guide_gff', cxt)
		if guide_gff is not None:
			st_args.extend(['-G', guide_gff])
			
		if self.get_parameter_value('stranded') is not None:
			st_args.append('--{}'.format(self.get_parameter_value_as_string('stranded')))
		
		if not self.get_parameter_value('cov_trim'):
			st_args.append('-t')
			
		if not self.get_parameter_value('mm_corr'):
			st_args.append('-u')
		
		if self.get_parameter_value('estimate'):
			st_args.append('-e')
			
		if self.get_parameter_value('verbose'):
			st_args.append('-v')
		
		scaffold_skip = self.get_parameter_value('skip_scaffold')	
		if scaffold_skip is not None and len(scaffold_skip) > 0:
			st_args.extend(['-x', ','.join(map(str, scaffold_skip))])	
		
		if self.get_parameter_value('out_ballgown'):
			st_args.append('-B')
			output_files.append(FileInfo(os.path.join(cxt.sample.dest, 'e_data.ctab'), FileContext.from_module_context(cxt, 'ballgown_exon_expression')))
			output_files.append(FileInfo(os.path.join(cxt.sample.dest, 'i_data.ctab'), FileContext.from_module_context(cxt, 'ballgown_intron_expression')))
			output_files.append(FileInfo(os.path.join(cxt.sample.dest, 't_data.ctab'), FileContext.from_module_context(cxt, 'ballgown_transcript_expression')))
			output_files.append(FileInfo(os.path.join(cxt.sample.dest, 'e2t.ctab'),	   FileContext.from_module_context(cxt, 'ballgown_exon_trans_map')))
			output_files.append(FileInfo(os.path.join(cxt.sample.dest, 'i2t.ctab'),    FileContext.from_module_context(cxt, 'ballgown_intron_trans_map')))
			
		if self.get_parameter_value('out_abundance'):
			out_abd = os.path.join(cxt.sample.dest, cxt.sample.name+'.gene_abund.tab')
			st_args.extend(['-A', out_abd])
			output_files.append(FileInfo(out_abd, FileContext.from_module_context(cxt, 'gene_abundance')))
			
		if self.get_parameter_value('out_full_cov'):
			out_full_cov = os.path.join(cxt.sample.dest, cxt.sample.name+'.full_cov_refs.gtf')
			st_args.extend(['-C', out_full_cov])
			output_files.append(FileInfo(out_full_cov, FileContext.from_module_context(cxt, 'full_coverage_ref_transcripts')))
			
		
		#add the input alignments that we are quantifying/assembling
		alignments = self.resolve_input('alignments', cxt)
		for a in alignments:
			st_args.append(a.fullpath)

		
		#OK, we now have all the arguments setup, lets actually run StringTie
		cxt.log.write("\t-> Assembling transcripts with StringTie......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(st_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(st_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		
		return output_files
	#end run()

#end class StringTieQuant





class StringTieMerge(StringTieBase):
	def __init__(self, **kwargs):
		super_args = dict(name='StringTieMerge', short_description='Merge Transcript assemblies with StringTie.')
		super_args.update(**kwargs)
		super(StringTieMerge, self).__init__(**super_args)
	#end __init__()

	def _declare_resolvers(self):
		self._name_resolver('assemblies')
		self._name_resolver('guide_gff')
	#end __declare_resolvers()

	
	def run(self, cxt):		
		gtf_out = os.path.join(cxt.sample.dest, cxt.sample.name+'.gtf'),
		
		st_args = [
			self.get_parameter_value_as_string('stringtie_path'),
			'--merge',
			'-p', str(self.processors), #not sure if this is actually supported, but lets give it a shot!
			'-l', cxt.sample.name,
			'-o', gtf_out,
			#'-G', #<guide_gff>, # @todo: need to implement this!
			
			#optional args
			'-m', self.get_parameter_value_as_string('min_len'),
			'-c', self.get_parameter_value_as_string('min_cov'),
			'-F', self.get_parameter_value_as_string('min_fpkm'),
			'-T', self.get_parameter_value_as_string('min_tpm'),
			'-f', self.get_parameter_value_as_string('min_iso'),
		]
		
		if self.get_parameter_value('keep_intron'):
			st_args.append('-i')
			
		guide_gff = self.resolve_input('guide_gff', cxt)
		if guide_gff is not None:
			st_args.extend(['-G', guide_gff])
		
		#add the input transcript assemblies that we are merging
		assemblies = self.resolve_input('assemblies', cxt)
		for a in assemblies:
			st_args.append(a.fullpath)
		
		#OK, we now have all the arguments setup, lets actually run StringTie
		cxt.log.write("\t-> Merging transcript assemblies with StringTie......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(st_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(st_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return [FileInfo(gtf_out, FileContext.from_module_context(cxt, 'merged_transcript_assembly'))]
		
	#end run()
#end class StringTieMerge
