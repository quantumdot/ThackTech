import os
import shutil 
import subprocess
from ThackTech.Pipelines import PipelineModule


class SortBam(PipelineModule):
	
	def __init__(self, **kwargs):
		super(SortBam, self).__init__('SortBam', 'Sort BAM file', **kwargs)
		self._name_resolver('alignments')
	#end __init__()

	def tool_versions(self):
		return {
			'samtools': subprocess.check_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		bam = self.resolve_input('alignments', cxt.sample)
		sorted_bam = os.path.splitext(bam)[0]+'_sorted'
		
		sort_cmd = ['samtools', 'sort', '-@', str(self.processors), bam, sorted_bam]
		
		#Sort our BAM file by genomic location: -@ => multithreading!
		cxt.log.write('-> Sorting BAM "%s" (using %d processor%s)...\n' % (bam, self.processors, ('s' if self.processors > 1 else '')))
		cxt.log.write("-> "+subprocess.check_output('samtools 2>&1 | grep Version', shell=True, stderr=subprocess.STDOUT)+"")
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(sort_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		
		self._run_subprocess(sort_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		#move the sorted BAM to overwrite the unsorted BAM
		os.remove(bam) #first remove the origional unsorted bam, otherwise the move will choke!
		shutil.move(sorted_bam+'.bam', bam) #mv "${dest}${alnname}_sorted.bam" "${dest}${alnname}.bam"
	#end run()
#end class GeneratePseudoreplicates