import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class IndexBam(PipelineModule):
	
	def __init__(self):
		super(IndexBam, self).__init__('IndexBam', 'Index BAM file')
		self._name_resolver('alignments')
	#end __init__()
	
	def tool_versions(self):
		return {
			'samtools': subprocess.check_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		bam = self.resolve_input('alignments', cxt.sample)
		
		index_cmd = ['samtools', 'index', bam]
		
		cxt.log.write('-> Indexing BAM "%s"...\n' % (bam,))
		cxt.log.write("-> "+subprocess.check_output('samtools 2>&1 | grep Version', shell=True, stderr=subprocess.STDOUT)+"")
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(index_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		
		self._run_subprocess(index_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
	#end run()
#end class GeneratePseudoreplicates