import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MergeBams(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MergeBams', short_description='Merge BAM files')
		super_args.update(**kwargs)
		super(MergeBams, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('samtools_path', str, 'samtools', desc="Path to samtools"))
		self.add_parameter(ModuleParameter('postfix', str, None, nullable=True, desc="Postfix to apply to merged bam filename"))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('alignments')
	#end __declare_resolvers()

	
	def tool_versions(self):
		return {
			'samtools': self._call_output(self.get_parameter_value('samtools_path')+" 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		outdir = os.path.join(cxt.sample.dest, 'merged_bam')
		filetools.ensure_dir(outdir)
		outbase = os.path.join(outdir, cxt.sample.name)
		if self.get_parameter_value('postfix') is not None:
			outbase += "."+self.get_parameter_value_as_string('postfix')
		outfile = outbase+'.merged.bam'
		
		samtools_cmd = [ self.get_parameter_value('samtools_path'), 'merge', outfile ] + self.resolve_input('alignments', cxt)
		
		cxt.log.write("-> Merging bam files with samtools......\n")
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(samtools_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(samtools_cmd, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return {
			'merged_bam': outfile
		}
	#end run()
#end class GeneratePseudoreplicates