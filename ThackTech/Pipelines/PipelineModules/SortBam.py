import os
import shutil 
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class SortBam(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='SortBam', short_description='Sort BAM file')
		super_args.update(**kwargs)
		super(SortBam, self).__init__(**super_args)
		
		self.add_parameter(ModuleParameter('sort_name', bool, False, desc="Sort by read name rather than coordinate."))
		self.add_parameter(ModuleParameter('overwrite', bool, True, desc="Overwrite the source, non-sorted, BAM."))
		
		self._name_resolver('alignments')
	#end __init__()

	def tool_versions(self):
		return {
			'samtools': self._call_output("samtools 2>&1 | perl -ne 'if(m/Version: ([\d\.-\w]+)/){ print $1; }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		bam = self.resolve_input('alignments', cxt)
		if bam is None:
			cxt.log.write("No bam provided, exiting module\n")
			return
		
		sorted_bam = os.path.join(cxt.sample.dest, bam.basename_with_ext('sorted'))
		
		sort_cmd = ['samtools', 'sort', '-@', str(self.processors), bam.fullpath, sorted_bam]
		
		if self.get_parameter_value('sort_name'):
			sort_cmd.insert(2, '-n')
		
		#Sort our BAM file by genomic location: -@ => multithreading!
		cxt.log.write('-> Sorting BAM "%s" (using %d processor%s)...\n' % (bam, self.processors, ('s' if self.processors > 1 else '')))
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(sort_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		
		self._run_subprocess(sort_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		if self.get_parameter_value('overwrite'):
			#move the sorted BAM to overwrite the unsorted BAM
			os.remove(bam) #first remove the origional unsorted bam, otherwise the move will choke!
			shutil.move(sorted_bam+'.bam', bam) #mv "${dest}${alnname}_sorted.bam" "${dest}${alnname}.bam"
		else:
			return [FileInfo(sorted_bam+'.bam', FileContext.from_module_context(cxt, "sorted_bam"))]
	#end run()
#end class GeneratePseudoreplicates