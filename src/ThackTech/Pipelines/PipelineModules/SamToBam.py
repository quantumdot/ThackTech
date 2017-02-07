import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class SamToBam(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'SamToBam', 'Convert SAM to BAM')
		self._name_resolver('sam')
	#end __init__()

	def supported_types(self):
		return ['sam']
	#end supported_types()
	
	def load_modules(self, logfile):
		pass
		# try:
			# #from ThackTech.Pipelines import lmodHelper
			# #logfile.write(lmodHelper.module("load", "samtools/0.1.19"))
			# from env_modules_python import module
			# logfile.write(module("load", "samtools/0.1.19"))
			# logfile.flush()
		# except:
			# pass
	#end load_modules()
	
	def run(self, sample, logfile):
		logfile.write("\t-> Postprocessing with samtools...\n")
		#self.load_modules(logfile)
		#print os.environ
		sam = self.resolve_input('sam', sample)
		bam = os.path.splitext(sam)[0]+'.bam'

		#convert SAM to BAM: -b => output BAM; -S => input is SAM; -@ => multithreading!
		logfile.write("\t-> Converting SAM to BAM (using %d processor%s)...\n" % (self.processors, ('s' if self.processors > 1 else '')))
		logfile.flush()
		self._run_subprocess(['samtools', 'view', '-b', '-S', '-@', str(self.processors), '-o', bam, sam])#samtools view -bS -@ $nump "${dest}${alnname}.sam" > "${dest}${alnname}.bam"
		
		#remove the SAM file as it is no longer needed
		os.remove(sam)
		
		#Sort our BAM file by genomic location: -@ => multithreading!
		logfile.write("\t-> Sorting BAM (using %d processor%s)...\n" % (self.processors, ('s' if self.processors > 1 else '')))
		logfile.flush()
		sorted_bam = os.path.splitext(bam)[0]+'_sorted'
		self._run_subprocess(['samtools', 'sort', '-@', str(self.processors), bam, sorted_bam])#samtools sort -@ $nump "${dest}${alnname}.bam" "${dest}${alnname}_sorted"
		
		#move the sorted BAM to overwrite the unsorted BAM
		os.remove(bam) #first remove the origional unsorted bam, otherwise the move will choke!
		self._run_subprocess(['mv', sorted_bam+'.bam', bam])
		#shutil.move() #mv "${dest}${alnname}_sorted.bam" "${dest}${alnname}.bam"
		
		#index our sorted BAM file
		logfile.write("\t-> Indexing BAM...\n")
		logfile.flush()
		self._run_subprocess(['samtools', 'index', bam])#samtools index "${dest}${alnname}.bam"
		
		return {
			'bam': bam,
			'bam_idx': bam+'.bai'
		}
	#end run()
#end class SamToBam