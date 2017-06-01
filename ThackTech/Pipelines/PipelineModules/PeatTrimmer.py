import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class PeatTrimmer(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='PEAT', short_description='Trim reads with PEAT')
		super_args.update(**kwargs)
		super(PeatTrimmer, self).__init__(**super_args)
	#end __init__()
	
	def __declare_parameters(self):
		self.add_parameter(ModuleParameter('l', int, 30, desc="Minimum gene fragment length, i.e. the fragment length for reverse complement check"))
		self.add_parameter(ModuleParameter('r', float, 0.3, desc="Mismatch rate applied in first stage reverse complement scan"))
		self.add_parameter(ModuleParameter('g', float, 0.6, desc="Mismatch rate applied in second stage gene portion check"))
		self.add_parameter(ModuleParameter('a', float, 0.4, desc="Mismatch rate applied in second stage adapter portion check"))
	#end __declare_parameters()
	
	def __declare_resolvers(self):
		self._name_resolver('fastq')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			#'peat': self._call_output(["PEAT", "-version"], stderr=subprocess.STDOUT) #peat doesnt seem to have a version option!
		}
	#end tool_versions()
	
	def run(self, cxt):
		if not cxt.sample.get_attribute('PE'):
			raise("PEAT only works with paired-ended data!")
		
		outdir = os.path.join(cxt.sample.dest, 'peat')
		filetools.ensure_dir(outdir) 

		#trimlog_loc = os.path.join(outdir, cxt.sample.name+'.trimlog')
		peat_args = [
			'PEAT', 'paired',
			'-n', str(self.processors),
			'-l', self.get_parameter_value_as_string('l'),
			'-r', self.get_parameter_value_as_string('r'),
			'-g', self.get_parameter_value_as_string('g'),
			'-a', self.get_parameter_value_as_string('a'),
			
		]
		out_files = []
		read_files = self.resolve_input('fastq', cxt)
		
		#Input Files
		peat_args.append('-1')
		peat_args.append([f for f in read_files if f.has_attribute_value("mate", 1)][0].fullpath)
		peat_args.append('-2')
		peat_args.append([f for f in read_files if f.has_attribute_value("mate", 2)][0].fullpath)
		#Output Files
			
		
		
		out_files.append(FileInfo(os.path.join(outdir, cxt.sample.name+'_R1.filtered.paired.fastq.gz'), 
						 FileContext.from_module_context(cxt, "filtered_paired_reads"),
						 mate=1, filtered=True, paired=True))
		out_files.append(FileInfo(os.path.join(outdir, cxt.sample.name+'_R1.filtered.unpaired.fastq.gz'), 
						 FileContext.from_module_context(cxt, "filtered_unpaired_reads"),
						 mate=1, filtered=True, paired=False))
		out_files.append(FileInfo(os.path.join(outdir, cxt.sample.name+'_R2.filtered.paired.fastq.gz'), 
						 FileContext.from_module_context(cxt, "filtered_paired_reads"),
						 mate=2, filtered=True, paired=True))
		out_files.append(FileInfo(os.path.join(outdir, cxt.sample.name+'_R2.filtered.unpaired.fastq.gz'), 
						 FileContext.from_module_context(cxt, "filtered_unpaired_reads"),
						 mate=2, filtered=True, paired=False))
		for f in out_files:
			peat_args.append(f.fullpath)
			

			
		
		
		cxt.log.write("\t-> Performing trimming of source sequences......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(peat_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(peat_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return out_files
	#end run()
#end class Trimmomatic	