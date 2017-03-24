import os
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class Trimmomatic(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'Trimmomatic', 'Trim reads with Trimmomatic')
		
		self.add_parameter(ModuleParameter('trim_adapt_fa_SE', str, '/mnt/ref/adapters/TruSeq3-SE.fa'))
		self.add_parameter(ModuleParameter('trim_adapt_fa_PE', str, '/mnt/ref/adapters/TruSeq3-PE-2.fa'))
		
		self.add_parameter(ModuleParameter('leading', int, 3, nullable=True, desc="Cut bases off the start of a read, if below a threshold quality"))
		self.add_parameter(ModuleParameter('trailing', int, 3, nullable=True, desc="Cut bases off the end of a read, if below a threshold quality"))
		self.add_parameter(ModuleParameter('sliding_window_width', int, 4))
		self.add_parameter(ModuleParameter('sliding_window_qthresh', int, 15))
		self.add_parameter(ModuleParameter('min_length', int, 25, nullable=True, desc="Drop the read if it is below a specified length"))
		
		#self.add_parameter(ModuleParameter('crop', int, None, nullable=True, desc="Cut the read to a specified length"))
		#self.add_parameter(ModuleParameter('head_crop', int, None, nullable=True, desc="Cut the specified number of bases from the start of the read"))
		
		self._name_resolver('fastq')
	#end __init__()
	
	def run(self, sample, logfile):
		outdir = os.path.join(sample.dest, 'trimmomatic')
		filetools.ensure_dir(outdir) 

		trimlog_loc = os.path.join(outdir, sample.name+'.trimlog')
		trimmomatic_args = [
			'Trimmomatic' + ('PE' if sample.get_attribute('PE') else 'SE'),
			'-threads', str(self.processors),
			'-trimlog', trimlog_loc
		]
		out_sequences = []
		read_files = self.resolve_input('fastq', sample)
		if sample.get_attribute('PE'):
			#Input Files
			trimmomatic_args.append(read_files[0])
			trimmomatic_args.append(read_files[1])
			#Output Files
			trimmomatic_args.append(os.path.join(outdir, sample.name+'_R1.filtered.paired.fastq.gz'))
			trimmomatic_args.append(os.path.join(outdir, sample.name+'_R1.filtered.unpaired.fastq.gz'))
			trimmomatic_args.append(os.path.join(outdir, sample.name+'_R2.filtered.paired.fastq.gz'))
			trimmomatic_args.append(os.path.join(outdir, sample.name+'_R2.filtered.unpaired.fastq.gz'))
			#files to return - properly paired and filtered sequences passing QC
			out_sequences.append(os.path.join(outdir, sample.name+'_R1.filtered.paired.fastq.gz'))
			out_sequences.append(os.path.join(outdir, sample.name+'_R2.filtered.paired.fastq.gz'))
		else:
			#Input Files
			trimmomatic_args.append(read_files[0])
			#Output Files
			trimmomatic_args.append(os.path.join(outdir, sample.name+'.filtered.fastq.gz'))
			#files to return - filtered sequences passing QC
			out_sequences.append(os.path.join(outdir, sample.name+'.filtered.fastq.gz'))
		#########################
		#  Trimming Parameters  #
		#########################
		if sample.get_attribute('PE'):
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
			
		
		
		logfile.write("\t-> Performing trimming of source sequences......")
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(trimmomatic_args))
		logfile.write("\n..............................................\n")
		logfile.flush()
		self._run_subprocess(trimmomatic_args)
		
		
		logfile.write("\t-> Compressing Trim Log......")
		logfile.flush()
		self._run_subprocess(['tar', 'cfz', trimlog_loc+'.tar.gz', trimlog_loc])
		
		
		return {
			'trimmed_reads': out_sequences,
			'trim_log':		 trimlog_loc
		}
	#end run()
#end class Trimmomatic	