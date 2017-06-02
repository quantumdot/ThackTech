import os
import subprocess
import re
import tempfile
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext
from ThackTech import filetools


class BowtieAlign(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='BowtieAlign', short_description='Alignment using Bowtie')
		super_args.update(**kwargs)
		super(BowtieAlign, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('bowtie_path', 		str, 	'bowtie',	desc="Path to the bowtie executable."))
		self.add_parameter(ModuleParameter('tryhard', 			bool,	False,	desc=""))
		self.add_parameter(ModuleParameter('best', 				bool, 	False,	desc="Report only the best possible singleton alignment."))
		self.add_parameter(ModuleParameter('unaligned', 		bool, 	False,	desc="Report reads that fail to align."))
		self.add_parameter(ModuleParameter('multimap', 			bool, 	False,	desc="Report reads that map to multiple locations in the reference."))
		self.add_parameter(ModuleParameter('max_align', 		int, 	None,	nullable=True, desc="If not none, maximum number of valid alignments to report."))
		self.add_parameter(ModuleParameter('chunkmbs', 			int, 	512,	desc="give more memory for searching.. prevents warnings and increases alignment rate especially for longer reads"))
		self.add_parameter(ModuleParameter('max_insert', 		int, 	1200,	desc="max insert size allowed for PE reads"))
		self.add_parameter(ModuleParameter('max_mismatches',	int, 	2,		desc="max allowed mismatches in the sead region, corresponds to -n"))
		self.add_parameter(ModuleParameter('pairtries', 		int, 	1000,	desc="number of tries for finding valid paired-end alignments"))
		self.add_parameter(ModuleParameter('additional_args', 	list, 	[],		desc="Additional arguments to pass to Bowtie"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastq')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'bowtie': self._call_output("bowtie --version 2>&1 | perl -ne 'if(m/.*bowtie.*version\s+([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing Bowtie....\n")
		cxt.log.flush()
		
		output_result = {}
		read_files = self.resolve_input('fastq', cxt)
				
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
			bowtiecmd += [
				'-1', [f for f in read_files if f.has_attribute_value("mate", 1)][0].fullpath, 
				'-2', [f for f in read_files if f.has_attribute_value("mate", 2)][0].fullpath
			]
		else:
			bowtiecmd.append(read_files[0].fullpath)
		#specify the destination SAM file
		output_result['sam'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.sam')
		bowtiecmd.append(output_result['sam'])

		
		#OK, we now have all the arguments setup, lets actually run bowtie
		cxt.log.write("\t-> Performing alignment with bowtie......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(bowtiecmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		tmpout = tempfile.NamedTemporaryFile()
		self._run_subprocess(bowtiecmd, stderr=subprocess.STDOUT, stdout=tmpout)
		tmpout.seek(0)
		bowtie_output = tmpout.read()
		tmpout.close()
		cxt.log.write(bowtie_output)
		output_result['align_stats'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.align_stats.tsv')
		self.parse_bowtie_output(cxt, bowtie_output, output_result['align_stats'])
		
		output_files = []
		for n, o in output_result.items():
			output_files.append(FileInfo(o, FileContext.from_module_context(cxt, n)))
		
		return output_files
	#end run()
	
	def parse_bowtie_output(self, cxt, logdata, destfilename):
		regex_items = [
			('Total Reads',		"# reads processed: (\d+)"),
			('Aligned Reads',	"# reads with at least one reported alignment: (\d+)"),
			('Unaligned Reads',	"# reads that failed to align: (\d+)"),
			('Percent Aligned Reads', lambda r: (float(r['Aligned Reads']) / float(r['Total Reads']))),
			('Percent Unaligned Reads', lambda r: (float(r['Unaligned Reads']) / float(r['Total Reads']))),
		]
		results = {}
		
		with open(destfilename, 'w') as destfile:
			destfile.write("Sample\t")
			destfile.write("\t".join([item[0] for item in regex_items]))
			destfile.write("\n")
			
			destfile.write("{}\t".format(cxt.sample.name))
			for item in regex_items:
				if isinstance(item[1], basestring):
					match = re.search(item[1], logdata, re.MULTILINE)
					if match is not None:
						results[item[0]] = int(match.group(1))
					else:
						results[item[0]] = 0
				else:
					results[item[0]] = item[1](results)
				
			destfile.write("\t".join(["{}".format(results[item[0]]) for item in regex_items]))
			destfile.write("\n")
	#end parse_bowtie_output
#end class BowtieAlign