import os
import re
import subprocess
import tempfile
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class Bowtie2Align(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Bowtie2Align', short_description='Alignment using Bowtie2')
		super_args.update(**kwargs)
		super(Bowtie2Align, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('bowtie2_path', 		str, 	'bowtie2',	desc="Path to the bowtie2 executable."))
		self.add_parameter(ModuleParameter('multimap', 			bool, 	False,	desc="Report reads that map to multiple locations in the reference."))
		self.add_parameter(ModuleParameter('max_insert', 		int, 	1200,	desc="--maxins: max insert size allowed for PE reads"))
		self.add_parameter(ModuleParameter('max_mismatches',	int, 	1,		desc="-N: max allowwd mismatches in the sead region"))
		self.add_parameter(ModuleParameter('unaligned', 		bool, 	False,	desc="--un: Report reads that fail to align."))
		self.add_parameter(ModuleParameter('no_unaligned', 		bool,	False,	desc="--no-unal: Suppress SAM records for reads that failed to align."))
		self.add_parameter(ModuleParameter('num_dva', 			int,	None, nullable=True, desc="Number of distinct, valid alignments to report. This corresponds to -k and -a. Leave None to disable, -1 for all, or a positive int."))
		self.add_parameter(ModuleParameter('additional_args', 	list, 	[],		desc="Additional arguments to pass to Bowtie2"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastq')
	#end __declare_resolvers()
	
	def load_modules(self):
		subprocess.call("module load bowtie2", shell=True)
	#end load_modules()
	
	def tool_versions(self):
		bt2p = self.get_parameter_value_as_string('bowtie2_path')
		return {
			'bowtie2': self._call_output(bt2p+" --version 2>&1 | perl -ne 'if(m/.*bowtie2.*version\s+([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		cxt.log.write("\t-> Preparing Bowtie2....\n")
		cxt.log.flush()
		
		output_result = {}
		read_files = self.resolve_input('fastq', cxt)
				
		#start constructing our bowtie2 arguments...
		bowtiecmd = [ 
			self.get_parameter_value_as_string('bowtie2_path'),
			'--threads', str(self.processors),
			'--time', 			#print time info			
			'-N', self.get_parameter_value_as_string('max_mismatches'), 			
		]
		
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += ['--maxins', 		self.get_parameter_value_as_string('max_insert')]
			
		
		if self.get_parameter_value('no_unaligned'):
			bowtiecmd.append('--no-unal')
			
		if self.get_parameter_value('num_dva') is not None:
			if self.get_parameter_value('num_dva') == -1:
				bowtiecmd.append('-a')
			else:
				bowtiecmd.extend(['-k', self.get_parameter_value_as_string('num_dva')])

		
		#check if we need to output reads that multimap
# 		if self.get_parameter_value('multimap'):
# 			if cxt.sample.get_attribute('PE'):
# 				output_result['multimap_1'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap_1.fastq')
# 				output_result['multimap_2'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap_2.fastq')
# 			else:
# 				output_result['multimap'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap.fastq')
# 			bowtiecmd += ['--max', os.path.join(cxt.sample.dest, cxt.sample.name+'_multimap.fastq')]
		
		#check if we need to output reads that fail to align
		if self.get_parameter_value('unaligned'):
			if cxt.sample.get_attribute('PE'):
				output_result['unaligned_1'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_1.fastq')
				output_result['unaligned_2'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned_2.fastq')
			else:
				output_result['unaligned'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')
			bowtiecmd += ['--un', os.path.join(cxt.sample.dest, cxt.sample.name+'_unaligned.fastq')]
		
		#hook in any additional args the user may want to supply to bowtie2
		bowtiecmd += self.get_parameter_value('additional_args')

			
		########################################
		#   All options must go BEFORE here!   #
		########################################
		
		#specify the reference genome index
		bowtiecmd += ['-x', cxt.sample.genome.get_index('Bowtie2Index')]
		#add the input file arguments
		if cxt.sample.get_attribute('PE'):
			bowtiecmd += [
				'-1', ",".join([f.fullpath for f in read_files if f.has_attribute_value("mate", 1)]), 
				'-2', ",".join([f.fullpath for f in read_files if f.has_attribute_value("mate", 2)])
			]
		else:
			bowtiecmd += [
				'-U', ",".join([f.fullpath for f in read_files])
			]
		#specify the destination SAM file
		output_result['sam'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.sam')
		bowtiecmd += ['-S', output_result['sam']]

		
		#OK, we now have all the arguments setup, lets actually run bowtie2
		cxt.log.write("\t-> Performing alignment with bowtie2......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join([str(p) for p in bowtiecmd]))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		tmpout = tempfile.NamedTemporaryFile()
		try:
			self._run_subprocess(bowtiecmd, stderr=subprocess.STDOUT, stdout=tmpout)
		finally:	
			tmpout.seek(0)
			bowtie_output = tmpout.read()
			tmpout.close()
			cxt.log.write(bowtie_output)
			cxt.log.flush()
			
		output_result['align_stats'] = os.path.join(cxt.sample.dest, cxt.sample.name+'.align_stats.tsv')
		self.parse_bowtie_output(cxt, bowtie_output, output_result['align_stats'])
		
		output_files = []
		for n, o in list(output_result.items()):
			output_files.append(FileInfo(o, FileContext.from_module_context(cxt, n)))
		
		return output_files
	#end run()
	
	def parse_bowtie_output(self, cxt, logdata, destfilename):
		if cxt.sample.get_attribute('PE'):
			regex_items = [
				('Count Total Reads',						  "(\d+) reads; of these:"),
				('Count Paired Reads', 						  "(\d+) \([\d\.%]+\) were paired; of these:"),
				('% Paired Reads',							  lambda r: (float(r['Count Total Reads']) / float(r['Count Paired Reads']))),
				('Count Paired Concordant Aligned 1 Time',	  "(\d+) \([\d\.%]+\) aligned concordantly exactly 1 time"),
				('% Paired Concordant Aligned 1 Time',		  lambda r: (float(r['Count Paired Concordant Aligned 1 Time']) / float(r['Count Paired Reads']))),
				('Count Paired Concordant Aligned > 1 Times', "(\d+) \([\d\.%]+\) aligned concordantly >1 times"),
				('% Paired Concordant Aligned > 1 Times', 	  lambda r: (float(r['Count Paired Concordant Aligned > 1 Times']) / float(r['Count Paired Reads']))),
				('% Overall Concordant Aligned', 			  lambda r: (float(r['Count Paired Concordant Aligned 1 Time'] + r['Count Paired Concordant Aligned > 1 Times']) / float(r['Count Paired Reads']))),
				('Count Paired Concordant Aligned 0 Times',   "(\d+) \([\d\.%]+\) aligned concordantly 0 times"),
				('Count Paired Discordant Aligned 1 Time', 	  "(\d+) \([\d\.%]+\) aligned discordantly 1 time"),
				('% Paired Discordant Aligned 1 Time', 		  lambda r: (float(r['Count Paired Discordant Aligned 1 Time']) / float(r['Count Paired Concordant Aligned 0 Times']))),
				('Count Paired Discordant Aligned > 1 Times', "(\d+) \([\d\.%]+\) aligned discordantly >1 times"),
				('% Paired Discordant Aligned > 1 Times', 	  lambda r: (float(r['Count Paired Discordant Aligned > 1 Times']) / float(r['Count Paired Concordant Aligned 0 Times']))),	
				('Count PE-Unaligned Pairs', 	  			  "(\d+) pairs aligned 0 times concordantly or discordantly; of these:"),
				('Count PE-Unaligned Reads', 	  			  "(\d+) mates make up the pairs; of these:"),
				('Count Unpaired Reads Aligned 0 Times', 	  "(\d+) \([\d\.%]+\) aligned 0 times"),
				('% Unpaired Reads Aligned 0 Times', 		  lambda r: (float(r['Count Unpaired Reads Aligned 1 Time']) / float(r['Count PE-Unaligned Reads']))),
				('Count Unpaired Reads Aligned 1 Time', 	  "(\d+) \([\d\.%]+\) aligned 0 times"),
				('% Unpaired Reads Aligned 1 Time', 		  lambda r: (float(r['Count Unpaired Reads Aligned 1 Time']) / float(r['Count PE-Unaligned Reads']))),
				('Count Unpaired Reads Aligned > 1 Times', 	  "(\d+) \([\d\.%]+\) aligned >1 times"),
				('% Unpaired Reads Aligned > 1 Times', 		  lambda r: (float(r['Count Unpaired Reads Aligned > 1 Time']) / float(r['Count PE-Unaligned Reads']))),
				
				('% Overall Alignment', lambda r: (( \
												 float((r['Count Paired Concordant Aligned 1 Time'] + r['Count Paired Concordant Aligned > 1 Times'] + r['Count Paired Discordant Aligned 1 Time'] + r['Count Paired Discordant Aligned > 1 Times']) * 2)\
												+ (r['Count Unpaired Reads Aligned 1 Time'] + r['Count Unpaired Reads Aligned > 1 Times'])) \
												/ float(r['Count Total Reads'] * 2) \
												))
			]
		else:
			regex_items = [
				('Total Reads',				"(\d+) reads; of these:"),
				('Reads Aligned 1 Time',	"# reads with at least one reported alignment: (\d+)"),
				('Reads Aligned >1 Times',	"# reads with at least one reported alignment: (\d+)"),
				('Reads Aligned 0 Times',	"# reads with at least one reported alignment: (\d+)"),
				('Percent Aligned Reads', 	lambda r: ((float(r['Reads Aligned 1 Time']) + float(r['Reads Aligned >1 Times'])) / float(r['Total Reads']))),
				('Percent Unaligned Reads', lambda r: (float(r['Reads Aligned 0 Times']) / float(r['Total Reads']))),
			]
		results = {}
		
		with open(destfilename, 'w') as destfile:
			destfile.write("Sample\t")
			destfile.write("\t".join([item[0] for item in regex_items]))
			destfile.write("\n")
			
			destfile.write("{}\t".format(cxt.sample.name))
			for item in regex_items:
				try:
					if isinstance(item[1], basestring):
						match = re.search(item[1], logdata, re.MULTILINE)
						if match is not None:
							results[item[0]] = int(match.group(1))
						else:
							results[item[0]] = 0
					else:
						results[item[0]] = item[1](results)
				except:
					results[item[0]] = 0
				
			destfile.write("\t".join(["{}".format(results[item[0]]) for item in regex_items]))
			destfile.write("\n")
	#end parse_bowtie_output
	
	
	
	
	
	
	
	
#end class BowtieAlign