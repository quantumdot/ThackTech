import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class HISAT2Align(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'HISAT2Align', 'Alignment using HISAT2')
		
		self.add_parameter(ModuleParameter('hisat2_path', 		str, 	'HISAT2',	desc="Path to the hisat2 executable."))
		
		self.add_parameter(ModuleParameter('unaligned', 		bool, 	False,	desc="Report reads that fail to align."))
		self.add_parameter(ModuleParameter('multimap', 			bool, 	False,	desc="Report reads that map to multiple locations in the reference."))
		self.add_parameter(ModuleParameter('max_align', 		int, 	None,	nullable=True, desc="If not none, maximum number of valid alignments to report."))
		self.add_parameter(ModuleParameter('max_insert', 		int, 	600,	desc="max insert size allowed for PE reads"))
		self.add_parameter(ModuleParameter('max_mismatches',	int, 	1,		desc="-N: max allowwd mismatches in the sead region"))
		self.add_parameter(ModuleParameter('no_unaligned', 		bool,	False,	desc="--no-unal: Suppress SAM records for reads that failed to align."))
		#self.add_parameter(ModuleParameter('tryhard', 			bool,	False,	desc=""))
		#self.add_parameter(ModuleParameter('best', 			bool, 	True,	desc="Report only the best possible singleton alignment."))
		#self.add_parameter(ModuleParameter('chunkmbs', 		int, 	512,	desc="give more memory for searching.. prevents warnings and increases alignment rate especially for longer reads"))
		#self.add_parameter(ModuleParameter('pairtries', 		int, 	1000,	desc="number of tries for finding valid paired-end alignments"))
		self.add_parameter(ModuleParameter('additional_args', 	list, 	[],		desc="Additional arguments to pass to Bowtie2"))

		self._name_resolver('fastq')
	#end __init__()

	def supported_types(self):
		return ['fastq']
	#end supported_types()
	
	def load_modules(self):
		subprocess.call("module load HISAT2", shell=True)
	#end load_modules()
	
	def show_version(self, handle=None, fancy=True):
		handle.write('\HISAT2 version is:\n------------------------------------------\n')
		handle.flush()
		proc = subprocess.Popen('HISAT2 --version', shell=True)
		proc.communicate()
		handle.write('------------------------------------------\n\n')
		handle.flush()
	#end show_version()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing HISAT2....\n")
		cxt.log.flush()
		
		output_result = {}
		read_files = self.resolve_input('fastq', cxt.sample)
				
		#start constructing our HISAT2 arguments...
		bowtiecmd = [ 
			self.get_parameter_value_as_string('hisat2_path'),
			'--threads', str(self.processors),
			'--time', 			#print time info			
			'-N', self.get_parameter_value_as_string('max_mismatches'), 			
		]
		
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += ['--maxins', 		self.get_parameter_value_as_string('max_insert')]
			
		if self.get_parameter_value('max_align') is not None:
			bowtiecmd += ['-m', self.get_parameter_value_as_string('max_align')]
		
		if self.get_parameter_value('no_unaligned'):
			bowtiecmd.append('--no-unal')
		#if self.get_parameter_value('best'):
		#	bowtiecmd.append('--best')
		#if self.get_parameter_value('tryhard'):
		#	bowtiecmd.append('--tryhard')
		
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
		
		#hook in any additional args the user may want to supply to HISAT2
		bowtiecmd += self.get_parameter_value('additional_args')

			
		########################################
		#   All options must go BEFORE here!   #
		########################################
		
		#specify the reference genome index
		bowtiecmd += ['-x', cxt.sample.genome.get_index('Bowtie2Index')]
		#add the input file arguments
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += ['-1', read_files[0], '-2', read_files[1]]
		else:
			bowtiecmd += ['-U', read_files[0]]
		#specify the destination SAM file
		output_result['sam'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.sam')
		bowtiecmd += ['-S', output_result['sam']]

		
		#OK, we now have all the arguments setup, lets actually run HISAT2
		cxt.log.write("\t-> Performing alignment with HISAT2......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bowtiecmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(bowtiecmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return output_result
	#end run()
#end class BowtieAlign