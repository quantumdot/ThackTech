import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class FastP(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='FastP', short_description='Read trimming and QC with FastP')
		super_args.update(**kwargs)
		super(FastP, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('fastp_path', str, 'fastp'))
		
		self.add_parameter(ModuleParameter('phred64', bool, False, desc="indicates the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)"))
		self.add_parameter(ModuleParameter('compression', int, 2, desc="compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest"))
		
		
		#Global Trimming
		self.add_parameter(ModuleParameter('trim_front1', int, 0, desc="trimming how many bases in front for read1"))
		self.add_parameter(ModuleParameter('trim_tail1', int, 0, desc="trimming how many bases in tail for read1"))
		self.add_parameter(ModuleParameter('trim_front2', int, None, nullable=True, desc="trimming how many bases in front for read2. If it's not specified, it will follow read1's settings"))
		self.add_parameter(ModuleParameter('trim_tail2', int, None, nullable=True, desc="trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings"))
		
		#polyG tail trimming
		self.add_parameter(ModuleParameter('polyG_trim', bool, None, nullable=True, desc="by default polyG tail trimming is automatically enabled for Illumina NextSeq/NovaSeq data. To force enable or disable, specify a value."))
		
		#per read cutting by quality score
		self.add_parameter(ModuleParameter('cut_by_quality5', bool, False, desc="enable per read cutting by quality in front (5')"))
		
		#UMI parameters
		self.add_parameter(ModuleParameter('enable_umi', bool, False, desc="enable unique molecular identifer (UMI) preprocessing"))
		self.add_parameter(ModuleParameter('umi_loc', str, None, nullable=True, choices=['index1','index2','read1','read2','per_index','per_read'], desc="specify the location of UMI"))
		self.add_parameter(ModuleParameter('umi_len', int, 0, desc="if the UMI is in read1/read2, its length should be provided"))
		self.add_parameter(ModuleParameter('umi_prefix', str, None, nullable=True, desc="if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG)."))
		
		
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
		
		out_files = [
			FileInfo(os.path.join(cxt.sample.dest, '{}.fastp.html'.format(cxt.sample.name)), FileContext.from_module_context(cxt, 'fastp_html_report')),
			FileInfo(os.path.join(cxt.sample.dest, '{}.fastp.json'.format(cxt.sample.name)), FileContext.from_module_context(cxt, 'fastp_json_report')),
		]

		fastp_args = [
			self.get_parameter_value_as_string('fastp_path'),
			'--thread', str(self.processors),
			'--html', out_files[0].fullpath,
			'--json', out_files[1].fullpath
		]
		
		read_files = self.resolve_input('fastq', cxt)
		if cxt.sample.get_attribute('PE'):
			#Input Files
			trimmomatic_args.append([f for f in read_files if f.has_attribute_value("mate", 1)][0].fullpath)
			trimmomatic_args.append([f for f in read_files if f.has_attribute_value("mate", 2)][0].fullpath)
			#Output Files

			
		
		
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