import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.ModuleParameter import ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class FastqScreen(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='FastqScreen', short_description='Screen FASTQ against reference genomes')
		super_args.update(**kwargs)
		super(FastqScreen, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('fastq_screen_path', str, 'fastq_screen', desc="Path to fastq_screen"))
		self.add_parameter(ModuleParameter('aligner', str, 'bowtie', choices=['bowtie', 'bowtie2', 'bwa'], desc="Aligner to use for the mapping"))
		self.add_parameter(ModuleParameter('conf', str, None, nullable=True, desc="Specify a location for the configuration (other than default)"))
		self.add_parameter(ModuleParameter('force', bool, True, desc="Do not terminate if output files already exist, instead overwrite the files."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('fastqs')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'fastq_screen': self._call_output(self.get_parameter_value('fastq_screen_path')+r" --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'fqscreen')
		filetools.ensure_dir(dest_dir)
		
		out_files = []
		fqscreen_args = [
			self.get_parameter_value('fastq_screen_path'),
			'--aligner', self.get_parameter_value_as_string('aligner'),
			'--threads', str(self.processors),
			'--outdir',  dest_dir
		]
		if self.get_parameter_value('conf') is not None:
			fqscreen_args += ['--conf', self.get_parameter_value_as_string('conf')]
		if self.get_parameter_value('force'):
			fqscreen_args.append('--force')
			
		for f in self.resolve_input('fastqs', cxt):
			fqscreen_args.append(f.fullpath)
			out_base = os.path.join(dest_dir, f.basename.replace(".fastq.gz", "_screen"))
			out_files.append(FileInfo("{}.html".format(out_base), FileContext.from_module_context(cxt, "FastqScreen_report")))
			out_files.append(FileInfo("{}.txt".format(out_base), FileContext.from_module_context(cxt, "FastqScreen_data")))
		
		

		cxt.log.write("Screening FASTQ for "+cxt.sample.name+" with Fastq Screen......\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(fqscreen_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(fqscreen_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return out_files
	#end run()
#end class BamFingerprint