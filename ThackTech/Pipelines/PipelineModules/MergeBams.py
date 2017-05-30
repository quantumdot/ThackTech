import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class MergeBams(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MergeBams', short_description='Merge BAM files')
		super_args.update(**kwargs)
		super(MergeBams, self).__init__(**super_args)
		
		self._name_resolver('alignments')
	#end __init__()

	
	def tool_versions(self):
		return {
			'samtools': self._call_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		outdir = os.path.join(cxt.sample.dest, 'merged_bam')
		filetools.ensure_dir(outdir)
		outfile = os.path.join(outdir, cxt.sample.name+'.merged.bam')
		
		samtools_cmd = [ 'samtools', 'merge', outfile ] + self.resolve_input('alignments', cxt)
		
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