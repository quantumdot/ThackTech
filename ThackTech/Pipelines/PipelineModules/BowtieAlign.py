import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class BowtieAlign(PipelineModule):
	
	def __init__(self, **kwargs):
		super(BowtieAlign, self).__init__('BowtieAlign', 'Alignment using Bowtie', **kwargs)
		
		self.add_parameter(ModuleParameter('bowtie_path', 		str, 	'bowtie',	desc="Path to the bowtie executable."))
		
		self.add_parameter(ModuleParameter('tryhard', 			bool,	False,	desc=""))
		self.add_parameter(ModuleParameter('best', 				bool, 	False,	desc="Report only the best possible singleton alignment."))
		self.add_parameter(ModuleParameter('unaligned', 		bool, 	False,	desc="Report reads that fail to align."))
		self.add_parameter(ModuleParameter('multimap', 			bool, 	False,	desc="Report reads that map to multiple locations in the reference."))
		self.add_parameter(ModuleParameter('max_align', 		int, 	None,	nullable=True, desc="If not none, maximum number of valid alignments to report."))
		self.add_parameter(ModuleParameter('chunkmbs', 			int, 	512,	desc="give more memory for searching.. prevents warnings and increases alignment rate especially for longer reads"))
		self.add_parameter(ModuleParameter('max_insert', 		int, 	600,	desc="max insert size allowed for PE reads"))
		self.add_parameter(ModuleParameter('max_mismatches',	int, 	2,		desc="max allowed mismatches in the sead region, corresponds to -n"))
		self.add_parameter(ModuleParameter('pairtries', 		int, 	1000,	desc="number of tries for finding valid paired-end alignments"))
		self.add_parameter(ModuleParameter('additional_args', 	list, 	[],		desc="Additional arguments to pass to Bowtie"))

		self._name_resolver('fastq')
	#end __init__()
	
	def tool_versions(self):
		return {
			'bowtie': subprocess.check_output("bowtie --version 2>&1 | perl -ne 'if(m/.*bowtie.*version\s+([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing Bowtie....\n")
		cxt.log.flush()
		
		output_result = {}
		read_files = self.resolve_input('fastq', cxt.sample)
				
		#start constructing our bowtie arguments...
		bowtiecmd = [ 
			self.get_parameter_value_as_string('bowtie_path'),
			'--threads', str(self.processors),
			'--sam', 			#output sam
			'--time', 			#print time info			
			'--chunkmbs', self.get_parameter_value_as_string('chunkmbs'),
			'-n', self.get_parameter_value_as_string('max_mismatches'), 			
		]
		
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += ['--pairtries', 	self.get_parameter_value_as_string('pairtries')]
			bowtiecmd += ['--maxins', 		self.get_parameter_value_as_string('max_insert')]
			
		if self.get_parameter_value('max_align') is not None:
			bowtiecmd += ['-m', self.get_parameter_value_as_string('max_align')]
		
		if self.get_parameter_value('best'):
			bowtiecmd.append('--best')
		if self.get_parameter_value('tryhard'):
			bowtiecmd.append('--tryhard')
		
		#check if we need to output reads that multimap
		if self.get_parameter_value('multimap'):
			if cxt.sample.get_attribute('PE'):
				output_result['multimap_1'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap_1.fastq')
				output_result['multimap_2'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap_2.fastq')
			else:
				output_result['multimap'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap.fastq')
			bowtiecmd += ['--max', os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap.fastq')]
		
		#check if we need to output reads that fail to align
		if self.get_parameter_value('unaligned'):
			if cxt.sample.get_attribute('PE'):
				output_result['unaligned_1'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_1.fastq')
				output_result['unaligned_2'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_2.fastq')
			else:
				output_result['unaligned'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')
			bowtiecmd += ['--un', os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')]
		
		#hook in any additional args the user may want to supply to bowtie
		bowtiecmd += self.get_parameter_value('additional_args')

			
		########################################
		#   All options must go BEFORE here!   #
		########################################
		
		#specify the reference genome index
		bowtiecmd.append(cxt.sample.genome.get_index('BowtieIndex'))
		#add the input file arguments
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += ['-1', read_files[0], '-2', read_files[1]]
		else:
			bowtiecmd.append(read_files[0])
		#specify the destination SAM file
		output_result['sam'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.sam')
		bowtiecmd.append(output_result['sam'])

		
		#OK, we now have all the arguments setup, lets actually run bowtie
		cxt.log.write("\t-> Performing alignment with bowtie......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bowtiecmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bowtiecmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return output_result
	#end run()
#end class BowtieAlign