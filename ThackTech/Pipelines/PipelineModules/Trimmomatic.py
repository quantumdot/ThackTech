import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class Trimmomatic(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='Trimmomatic', short_description='Trim reads with Trimmomatic')
		super_args.update(**kwargs)
		super(Trimmomatic, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('trim_adapt_fa_SE', str, '/mnt/ref/adapters/TruSeq3-SE.fa'))
		self.add_parameter(ModuleParameter('trim_adapt_fa_PE', str, '/mnt/ref/adapters/TruSeq3-PE-2.fa'))
		self.add_parameter(ModuleParameter("compress_trimlog", bool, True, desc="Determines if the trim log is compressed."))
		self.add_parameter(ModuleParameter('leading', int, 3, nullable=True, desc="Cut bases off the start of a read, if below a threshold quality"))
		self.add_parameter(ModuleParameter('trailing', int, 3, nullable=True, desc="Cut bases off the end of a read, if below a threshold quality"))
		self.add_parameter(ModuleParameter('sliding_window_width', int, 4))
		self.add_parameter(ModuleParameter('sliding_window_qthresh', int, 15))
		self.add_parameter(ModuleParameter('min_length', int, 25, nullable=True, desc="Drop the read if it is below a specified length"))
		#self.add_parameter(ModuleParameter('crop', int, None, nullable=True, desc="Cut the read to a specified length"))
		#self.add_parameter(ModuleParameter('head_crop', int, None, nullable=True, desc="Cut the specified number of bases from the start of the read"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastq')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'Trimmomatic': self._call_output(["TrimmomaticSE", "-version"], stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		outdir = os.path.join(cxt.sample.dest, 'trimmomatic')
		filetools.ensure_dir(outdir) 

		trimlog_loc = os.path.join(outdir, cxt.sample.name+'.trimlog')
		trimmomatic_args = [
			'Trimmomatic' + ('PE' if cxt.sample.get_attribute('PE') else 'SE'),
			'-threads', str(self.processors),
			'-trimlog', trimlog_loc
		]
		out_files = []
		read_files = self.resolve_input('fastq', cxt)
		if cxt.sample.get_attribute('PE'):
			#Input Files
			trimmomatic_args.append([f for f in read_files if f.has_attribute_value("mate", 1)][0].fullpath)
			trimmomatic_args.append([f for f in read_files if f.has_attribute_value("mate", 2)][0].fullpath)
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
				trimmomatic_args.append(f.fullpath)
			
# 			trimmomatic_args.append(os.path.join(outdir, cxt.sample.name+'_R1.filtered.paired.fastq.gz'))
# 			trimmomatic_args.append(os.path.join(outdir, cxt.sample.name+'_R1.filtered.unpaired.fastq.gz'))
# 			trimmomatic_args.append(os.path.join(outdir, cxt.sample.name+'_R2.filtered.paired.fastq.gz'))
# 			trimmomatic_args.append(os.path.join(outdir, cxt.sample.name+'_R2.filtered.unpaired.fastq.gz'))
# 			#files to return - properly paired and filtered sequences passing QC
# 			out_sequences.append(os.path.join(outdir, cxt.sample.name+'_R1.filtered.paired.fastq.gz'))
# 			out_sequences.append(os.path.join(outdir, cxt.sample.name+'_R2.filtered.paired.fastq.gz'))
		else:
			#Input Files
			trimmomatic_args.append(read_files[0].fullpath)
			#Output Files
			out_files.append(FileInfo(os.path.join(outdir, cxt.sample.name+'.filtered.fastq.gz'), 
							 FileContext.from_module_context(cxt, "filtered_reads"),
							 filtered=True))
			for f in out_files:
				trimmomatic_args.append(f.fullpath)
			#trimmomatic_args.append(os.path.join(outdir, cxt.sample.name+'.filtered.fastq.gz'))
			#files to return - filtered sequences passing QC
			#out_sequences.append(os.path.join(outdir, cxt.sample.name+'.filtered.fastq.gz'))
		#########################
		#  Trimming Parameters  #
		#########################
		if cxt.sample.get_attribute('PE'):
			#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clipthreshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
			trimmomatic_args.append('ILLUMINACLIP:%s:%d:%d:%d:%d:%s' % (self.get_parameter_value('trim_adapt_fa_PE'), 2, 30, 10, 2, 'true'))
		else:
			#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clipthreshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
			trimmomatic_args.append('ILLUMINACLIP:%s:%d:%d:%d:%d:%s' % (self.get_parameter_value('trim_adapt_fa_SE'), 2, 30, 10, 2, 'true'))
		
		if self.get_parameter_value('leading') is not None:
			#LEADING:<quality>
			trimmomatic_args.append('LEADING:%d' % (self.get_parameter_value('leading'),))
			
		if self.get_parameter_value('trailing') is not None:
			#TRAILING:<quality>
			trimmomatic_args.append('TRAILING:%d' % (self.get_parameter_value('trailing'),))
		#SLIDINGWINDOW:<windowSize>:<requiredQuality>
		trimmomatic_args.append('SLIDINGWINDOW:%d:%d' % (4,15))
		
		if self.get_parameter_value('min_length') is not None:
			#MINLEN:<length> 
			trimmomatic_args.append('MINLEN:%d' % (self.get_parameter_value('min_length'),))
			
		
		
		cxt.log.write("\t-> Performing trimming of source sequences......")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(trimmomatic_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(trimmomatic_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		if self.get_parameter_value("compress_trimlog"):
			cxt.log.write("\t-> Compressing Trim Log......\n")
			cxt.log.flush()
			self._run_subprocess(['tar', 'cfz', trimlog_loc+'.tar.gz', trimlog_loc])
			trimlog_loc = trimlog_loc+'.tar.gz'
		
		out_files.append(FileInfo(trimlog_loc, FileContext.from_module_context(cxt, "trimlog")))
		
		return out_files
	#end run()
#end class Trimmomatic	