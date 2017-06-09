import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter, FileInfo, FileContext


class IndexBam(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='IndexBam', short_description='Index BAM file')
		super_args.update(**kwargs)
		super(IndexBam, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('samtools_path', str, 'samtools', desc="Path to samtools"))
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
		bam = self.resolve_input('alignments', cxt)
		
		index_cmd = [self.get_parameter_value('samtools_path'), 'index', bam.fullpath]
		
		cxt.log.write('-> Indexing BAM "%s"...\n' % (bam.fullpath,))
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(index_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		
		self._run_subprocess(index_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		bam.companions.append(FileInfo(bam.fullpath+'.bai', FileContext.from_module_context(cxt, "bam_index")))
	#end run()
#end class GeneratePseudoreplicates