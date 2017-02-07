import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class IndexBam(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'IndexBam', 'Index BAM file')
		self._name_resolver('alignments')
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	def run(self, sample, logfile):
		bam = self.resolve_input('alignments', sample)
		
		index_cmd = ['samtools', 'index', bam]
		
		logfile.write('-> Indexing BAM "%s"...\n' % (bam,))
		logfile.write("-> "+subprocess.check_output('samtools 2>&1 | grep Version', shell=True, stderr=subprocess.STDOUT)+"")
		logfile.write("..............................................\n")
		logfile.write(" ".join(index_cmd))
		logfile.write("\n..............................................\n")
		logfile.flush()
		
		self._run_subprocess(index_cmd, stderr=subprocess.STDOUT, stdout=logfile)
	#end run()
#end class GeneratePseudoreplicates