import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MergeBams(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'MergeBams', 'Merge BAM files')
		self._name_resolver('alignments')
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	def tool_versions(self):
		return {
			'samtools': subprocess.check_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		outdir = os.path.join(cxt.sample.dest, 'merged_bam')
		filetools.ensure_dir(outdir)
		outfile = os.path.join(outdir, cxt.sample.name+'.merged.bam')
		
		samtools_cmd = [ 'samtools', 'merge', outfile ] + self.resolve_input('alignments', cxt.sample)
		
		cxt.log.write("-> Merging bam files with samtools......\n")
		cxt.log.write("-> "+subprocess.check_output('samtools 2>&1 | grep Version', shell=True, stderr=subprocess.STDOUT)+"")
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