import os
import shutil 
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class SortBam(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'SortBam', 'Sort BAM file')
		self._name_resolver('alignments')
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	def run(self, sample, logfile):
		bam = self.resolve_input('alignments', sample)
		sorted_bam = os.path.splitext(bam)[0]+'_sorted'
		
		sort_cmd = ['samtools', 'sort', '-@', str(self.processors), bam, sorted_bam]
		
		#Sort our BAM file by genomic location: -@ => multithreading!
		logfile.write('-> Sorting BAM "%s" (using %d processor%s)...\n' % (bam, self.processors, ('s' if self.processors > 1 else '')))
		logfile.write("-> "+subprocess.check_output('samtools 2>&1 | grep Version', shell=True, stderr=subprocess.STDOUT)+"")
		logfile.write("..............................................\n")
		logfile.write(" ".join(sort_cmd))
		logfile.write("\n..............................................\n")
		logfile.flush()
		
		self._run_subprocess(sort_cmd, stderr=subprocess.STDOUT, stdout=logfile)
		
		#move the sorted BAM to overwrite the unsorted BAM
		os.remove(bam) #first remove the origional unsorted bam, otherwise the move will choke!
		shutil.move(sorted_bam+'.bam', bam) #mv "${dest}${alnname}_sorted.bam" "${dest}${alnname}.bam"
	#end run()
#end class GeneratePseudoreplicates