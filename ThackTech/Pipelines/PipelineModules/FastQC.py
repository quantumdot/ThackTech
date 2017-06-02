import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.ModuleParameter import ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class FastQC(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='FastQC', short_description='FASTQ Quality Metrics')
		super_args.update(**kwargs)
		super(FastQC, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('kmers', int, 7, desc="Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10."))
		self.add_parameter(ModuleParameter('limits', str, None, nullable=True, desc="Specify a location for the limits file (--limits option in fastqc)"))
		self.add_parameter(ModuleParameter('adapters', str, None, nullable=True, desc="Specify a location for the adapters file (--adapters option in fastqc)"))
		self.add_parameter(ModuleParameter('contaminants', str, None, nullable=True, desc="Specify a location for the contaminants file (--contaminants option in fastqc)"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastqs')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'fastqc': self._call_output("fastqc --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'fastqc')
		filetools.ensure_dir(dest_dir)
		
		out_files = []
		fastqc_args = [
			'fastqc',
			'--quiet', #disable progress reporting
			'--threads', str(self.processors),
			'--outdir',  dest_dir,
			'--kmers', self.get_parameter_value_as_string('kmers')
		]
		if self.get_parameter_value('limits') is not None:
			fastqc_args += ['--limits', self.get_parameter_value_as_string('limits')]
			
		if self.get_parameter_value('adapters') is not None:
			fastqc_args += ['--adapters', self.get_parameter_value_as_string('adapters')]
		
		if self.get_parameter_value('contaminants') is not None:
			fastqc_args += ['--contaminants', self.get_parameter_value_as_string('contaminants')]
			
			
		for f in self.resolve_input('fastqs', cxt):
			fastqc_args.append(f.fullpath)
			of = FileInfo(os.path.join(dest_dir, "{}_fastqc.html".format(f.filename)), FileContext.from_module_context(cxt, "FastQC_report"))
			of.companions.append(FileInfo(os.path.join(dest_dir, "{}_fastqc.zip".format(f.filename)), FileContext.from_module_context(cxt, "FastQC_report_support_files")))
			out_files.append(of)
		

		cxt.log.write("Running FastQC for "+cxt.sample.name+"......\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(fastqc_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(fastqc_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return out_files
	#end run()
#end class BamFingerprint