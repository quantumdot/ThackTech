import os
import subprocess
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class SamToBam(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='SamToBam', short_description='Convert SAM to BAM')
		super_args.update(**kwargs)
		super(SamToBam, self).__init__(**super_args)
		
		self._name_resolver('sam')
	#end __init__()
	
	def tool_versions(self):
		return {
			'samtools': self._call_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def load_modules(self, cxt):
		pass
		# try:
			# #from ThackTech.Pipelines import lmodHelper
			# #cxt.log.write(lmodHelper.module("load", "samtools/0.1.19"))
			# from env_modules_python import module
			# cxt.log.write(module("load", "samtools/0.1.19"))
			# cxt.log.flush()
		# except:
			# pass
	#end load_modules()
	
	def run(self, cxt):
		cxt.log.write("\t-> Postprocessing with samtools...\n")
		#self.load_modules(cxt.log)
		#print os.environ
		sam = self.resolve_input('sam', cxt)
		bam = os.path.splitext(sam.fullpath)[0]+'.bam'

		#convert SAM to BAM: -b => output BAM; -S => input is SAM; -@ => multithreading!
		cxt.log.write("\t-> Converting SAM to BAM (using %d processor%s)...\n" % (self.processors, ('s' if self.processors > 1 else '')))
		cxt.log.flush()
		self._run_subprocess(['samtools', 'view', '-b', '-S', '-@', str(self.processors), '-o', bam, sam.fullpath])#samtools view -bS -@ $nump "${dest}${alnname}.sam" > "${dest}${alnname}.bam"
		
		#remove the SAM file as it is no longer needed
		cxt.sample.remove_file(sam)
		os.remove(sam.fullpath)
		
		#Sort our BAM file by genomic location: -@ => multithreading!
		cxt.log.write("\t-> Sorting BAM (using %d processor%s)...\n" % (self.processors, ('s' if self.processors > 1 else '')))
		cxt.log.flush()
		sorted_bam = os.path.splitext(bam)[0]+'_sorted'
		self._run_subprocess(['samtools', 'sort', '-@', str(self.processors), bam, sorted_bam])#samtools sort -@ $nump "${dest}${alnname}.bam" "${dest}${alnname}_sorted"
		
		#move the sorted BAM to overwrite the unsorted BAM
		os.remove(bam) #first remove the origional unsorted bam, otherwise the move will choke!
		self._run_subprocess(['mv', sorted_bam+'.bam', bam])
		#shutil.move() #mv "${dest}${alnname}_sorted.bam" "${dest}${alnname}.bam"
		
		#index our sorted BAM file
		cxt.log.write("\t-> Indexing BAM...\n")
		cxt.log.flush()
		self._run_subprocess(['samtools', 'index', bam])#samtools index "${dest}${alnname}.bam"
		
		
		outbam = FileInfo(bam, FileContext.from_module_context(cxt, "bam"))
		outbam.companions.append(FileInfo(bam+'.bai', FileContext.from_module_context(cxt, "bam_index")))
		return [outbam]
	#end run()
#end class SamToBam