import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class HISAT2Align(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='HISAT2Align', short_description='Alignment using HISAT2')
		super_args.update(**kwargs)
		super(HISAT2Align, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('hisat2_path', 		str, 	'HISAT2',	desc="Path to the hisat2 executable."))
		
		#Input:
		self.add_parameter(ModuleParameter('skip', int, None, nullable=True, desc="skip the first <int> reads/pairs in the input"))
		self.add_parameter(ModuleParameter('upto', int, None, nullable=True, desc="stop after first <int> reads/pairs"))
		
		#Alignment:
		self.add_parameter(ModuleParameter('n_ceil', str, 'L,0,0.15', desc="func for max # non-A/C/G/Ts permitted in alignment"))
		self.add_parameter(ModuleParameter('ignore_quals', bool, False, desc="treat all quality values as 30 on Phred scale"))
		self.add_parameter(ModuleParameter('nofw', bool, False, desc="do not align forward (original) version of read"))
		self.add_parameter(ModuleParameter('norc', bool, False, desc="do not align reverse-complement version of read"))
		
		#Spliced Alignment:
		self.add_parameter(ModuleParameter('no_spliced_alignment', bool, False, desc="disable spliced alignment"))
		self.add_parameter(ModuleParameter('pen_cansplice', int, 0, desc="penalty for a canonical splice site"))
		self.add_parameter(ModuleParameter('pen_noncansplice', int, 12, desc="penalty for a non_canonical splice site"))
		self.add_parameter(ModuleParameter('pen_canintronlen', str, 'G,-8,1', desc="penalty for long introns with canonical splice sites"))
		self.add_parameter(ModuleParameter('pen_noncanintronlen', str, 'G,-8,1', desc="penalty for long introns with noncanonical splice sites"))
		self.add_parameter(ModuleParameter('min_intronlen', int, 20, desc="minimum intron length"))
		self.add_parameter(ModuleParameter('max_intronlen', int, 500000, desc="maximum intron length"))
		self.add_parameter(ModuleParameter('no_temp_splicesite', bool, False, desc="disable the use of splice sites found"))
		self.add_parameter(ModuleParameter('rna_strandness', str, 'unstranded', desc="specify strand-specific information"))
		self.add_parameter(ModuleParameter('tmo', bool, False, desc="reports only those alignments within known transcriptome"))
		self.add_parameter(ModuleParameter('dta', str, None, choices=['stringtie', 'cufflinks'], nullable=True, desc="Specify which, if any, downstream, transcript assemblers output should be tailored to. see --dta* option in HISAT2"))
		self.add_parameter(ModuleParameter('avoid_pseudogene', bool, False, desc="tries to avoid aligning reads to pseudogenes; experimental!"))
		self.add_parameter(ModuleParameter('no_templatelen_adjustment', bool, False, desc="disables template length adjustment for RNA_seq reads"))
		self.add_parameter(ModuleParameter('write_novel_splicesites', bool, False, desc="report a list of splice sites"))
		
		#Scoring
		self.add_parameter(ModuleParameter('pen_mismatch_min', int, 2, desc="min penalty for mismatch; lower qual = lower penalty"))
		self.add_parameter(ModuleParameter('pen_mismatch_max', int, 6, desc="max penalty for mismatch; lower qual = lower penalty"))
		self.add_parameter(ModuleParameter('pen_sclip_min', int, 1, desc="min penalty for soft-clipping; lower qual = lower penalty"))
		self.add_parameter(ModuleParameter('pen_sclip_max', int, 2, desc="max penalty for soft-clipping; lower qual = lower penalty"))
		self.add_parameter(ModuleParameter('no_softclip', bool, False, desc="no soft-clipping"))
		self.add_parameter(ModuleParameter('pen_n', int, 1, desc="penalty for non-A/C/G/Ts in read/ref"))
		self.add_parameter(ModuleParameter('pen_read_gap_open', int, 5, desc="read gap open penalty"))
		self.add_parameter(ModuleParameter('pen_read_gap_extn', int, 3, desc="read gap extend penalty"))
		self.add_parameter(ModuleParameter('pen_ref_gap_open', int, 5, desc="reference gap open penalty"))
		self.add_parameter(ModuleParameter('pen_ref_gap_extn', int, 3, desc="reference gap extend penalty"))
		self.add_parameter(ModuleParameter('score_min', str, 'L,0.0,-0.2', desc="min acceptable alignment score with respect to read length"))
		
		#Reporting
		#-k
		
		#Paired-end:
		self.add_parameter(ModuleParameter('minins', int, 0, desc="minimum fragment length, only valid with --no-spliced-alignment"))
		self.add_parameter(ModuleParameter('maxins', int, 500, desc="maximum fragment length, only valid with --no-spliced-alignment"))
		self.add_parameter(ModuleParameter('mate_orient', str, 'fr', choices=['fr', 'rf', 'ff'], desc='-1, -2 mates align fw/rev, rev/fw, fw/fw'))
		self.add_parameter(ModuleParameter('no_mixed', bool, False, desc="suppress unpaired alignments for paired reads"))
		self.add_parameter(ModuleParameter('no_discordant', bool, False, desc="suppress discordant alignments for paired reads"))
		
		#Performance
		self.add_parameter(ModuleParameter('reorder', bool, False, desc="force SAM output order to match order of input reads"))
		self.add_parameter(ModuleParameter('offrate', int, None, nullable=True, desc="If not None, override offrate of index; must be >= index's offrate"))
		
		#Other
		self.add_parameter(ModuleParameter('unaligned', 		bool, 	False,	desc="Report reads that fail to align."))
		self.add_parameter(ModuleParameter('additional_args', 	list, 	[],		desc="Additional arguments to pass to Bowtie2"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastq')
		self._name_resolver('known_splicesites')
		self._name_resolver('novel_splicesites')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'HISAT2': self._call_output("hisat2 --version 2>&1 | perl -ne 'if(m/.*hisat2.*version\s+([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing HISAT2....\n")
		cxt.log.flush()
		
		output_result = {}
		read_files = self.resolve_input('fastq', cxt)
		output_result['sam'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.sam')
				
		#start constructing our HISAT2 arguments...
		bowtiecmd = [ 
			self.get_parameter_value_as_string('hisat2_path'),
			'-x', cxt.sample.genome.get_index('Hisat2Index'),
		]
		if cxt.sample.get_attribute('PE'):
			bowtiecmd.extend([
				'-1', ",".join([f.fullpath for f in read_files if f.has_attribute_value("mate", 1)]), 
				'-2', ",".join([f.fullpath for f in read_files if f.has_attribute_value("mate", 2)])
			])
		else:
			bowtiecmd.extend(['-U', ",".join([f.fullpath for f in read_files])])
		
		bowtiecmd.extend([
			'-S', output_result['sam'],
			'--threads', str(self.processors),
			'--time', 			#print time info			
			'--n-ceil', self.get_parameter_value_as_string('n_ceil'),
			
			#Scoring:
			'--mp', "{max},{min}".format(min=self.get_parameter_value('pen_mismatch_min'), max=self.get_parameter_value('pen_mismatch_min')),
  			'--sp', "{max},{min}".format(min=self.get_parameter_value('pen_sclip_min'), max=self.get_parameter_value('pen_sclip_max')),
    		
			'--np', self.get_parameter_value('pen_n'),
			'--rdg', "{opn},{ext}".format(opn=self.get_parameter_value('pen_read_gap_open'), ext=self.get_parameter_value('pen_read_gap_extn')),
			'--rfg', "{opn},{ext}".format(opn=self.get_parameter_value('pen_ref_gap_open'), ext=self.get_parameter_value('pen_ref_gap_extn')),
			'--score-min',  self.get_parameter_value_as_string('score_min'),
		])

		
		if self.get_parameter_value('skip') is not None:
			bowtiecmd.extend(['--skip', self.get_parameter_value_as_string('skip')])
			
		if self.get_parameter_value('upto') is not None:
			bowtiecmd.extend(['--upto', self.get_parameter_value_as_string('upto')])
		
		if self.get_parameter_value('reorder'):
			bowtiecmd.append('--reorder')
			
		if self.get_parameter_value('offrate') is not None:
			bowtiecmd.extend(['--offrate', self.get_parameter_value_as_string('offrate')])
		
		if self.get_parameter_value('no_softclip'):
			bowtiecmd.append('--no-softclip')
		
		if self.get_parameter_value('ignore_quals'):
			bowtiecmd.append('--ignore-quals')
			
		if self.get_parameter_value('nofw'):
			bowtiecmd.append('--nofw')
			
		if self.get_parameter_value('nofw'):
			bowtiecmd.append('--norc')
		
		
		if self.get_parameter_value('no_spliced_alignment'):
			bowtiecmd.append('--no-spliced-alignment')
		else:
			#spliced alignment enabled
			bowtiecmd.extend([
				'--pen-cansplice', self.get_parameter_value_as_string('pen_cansplice'),
				'--pen-noncansplice', self.get_parameter_value_as_string('pen_noncansplice'),
				'--pen-canintronlen', self.get_parameter_value_as_string('pen_canintronlen'),
				'--pen-noncanintronlen', self.get_parameter_value_as_string('pen_noncanintronlen'),
				'--min-intronlen', self.get_parameter_value_as_string('min_intronlen'),
				'--max-intronlen', self.get_parameter_value_as_string('max_intronlen'),
				'--rna-strandness', self.get_parameter_value_as_string('rna_strandness'),
			])
			if self.get_parameter_value('no_temp_splicesite'):
				bowtiecmd.append('--no-temp-splicesite')
				
			if self.get_parameter_value('tmo'):
				bowtiecmd.append('--tmo')
				
			if self.get_parameter_value('avoid_pseudogene'):
				bowtiecmd.append('--avoid-pseudogene')
				
			if self.get_parameter_value('no_templatelen_adjustment'):
				bowtiecmd.append('--no-templatelen-adjustment')
			
			dta_val = self.get_parameter_value('dta')
			if dta_val is not None:
				if dta_val == 'cufflinks':
					bowtiecmd.append('--dta-cufflinks')
				else:
					bowtiecmd.append('--dta')
			
			if self.get_parameter_value('write_novel_splicesites'):
				nss = FileInfo(os.path.join(cxt.sample.dest, 'novel_splicesites.txt'), FileContext.from_module_context(cxt, 'novel_splicesites'))
				bowtiecmd.extend(['--novel-splicesite-infile', nss.fullpath])
				
			# @todo: implement --novel-splicesite-infile <path> and --novel-splicesite-infile <path>
			
		#end if no_spliced_alignment
		
		
		if cxt.sample.get_attribute('PE'):
			bowtiecmd.extend([
				'--minins', self.get_parameter_value_as_string('minins'),
				'--maxins', self.get_parameter_value_as_string('maxins'),
				'--{}'.format(self.get_parameter_value_as_string('mate_orient'))
			])
			if self.get_parameter_value('no_mixed'):
				bowtiecmd.append('--no-mixed')
				
			if self.get_parameter_value('no_discordant'):
				bowtiecmd.append('--no-discordant') 
		#end PE options	
		
		#check if we need to output reads that fail to align
		if self.get_parameter_value('unaligned'):
			if cxt.sample.get_attribute('PE'):
				output_result['unaligned_1'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_1.fastq')
				output_result['unaligned_2'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_2.fastq')
			else:
				output_result['unaligned'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')
			bowtiecmd.extend(['--un', os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')])
		
		#hook in any additional args the user may want to supply to HISAT2
		bowtiecmd.extend(self.get_parameter_value('additional_args'))


		
		#OK, we now have all the arguments setup, lets actually run HISAT2
		cxt.log.write("\t-> Performing alignment with HISAT2......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(str(bowtiecmd))
		cxt.log.write(" ".join(bowtiecmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bowtiecmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return output_result
	#end run()
#end class HISAT2Align