import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.ModuleParameter import ModuleParameter
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class BamFingerprint(PipelineModule):
	
	def __init__(self, **kwargs):
		super(BamFingerprint, self).__init__('BamFingerprint', 'Bam Fingerprint', **kwargs)
		
		self.add_parameter(ModuleParameter('plotFileFormat', str, 'pdf', choices=['png', 'pdf', 'svg', 'eps'], desc="Plot output format"))
		self.add_parameter(ModuleParameter('outputRawCounts', bool, True, desc="Output raw count data"))
		
		self._name_resolver('bams')
	#end __init__()
	
	def tool_versions(self):
		return {
			'plotFingerprint': self._call_output("plotFingerprint --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'fingerprint')
		filetools.ensure_dir(dest_dir)
		outformat = self.get_parameter_value_as_string('plotFileFormat')
		out_files = []
		fingerprint_args = [
			#'bamFingerprint',
			'plotFingerprint',	#renamed in deepTools 2.0
			'--numberOfProcessors', str(self.processors),
			'--plotFile',  os.path.join(dest_dir, '{}.fingerprint.{}'.format(cxt.sample.name, outformat)),
			'--plotFileFormat', outformat
		]
		out_files.append(FileInfo( os.path.join(dest_dir, '{}.fingerprint.{}'.format(cxt.sample.name, outformat)),
						 FileContext.from_module_context(cxt, 'figure')))
		
		if self.get_parameter_value('outputRawCounts'):
			fingerprint_args += [
				'--outRawCounts', 
				os.path.join(dest_dir, '{}.fingerprint.rawcounts.txt'.format(cxt.sample.name))
			]
			out_files.append(FileInfo(os.path.join(dest_dir, '{}.fingerprint.rawcounts.txt'.format(cxt.sample.name)),
							 FileContext.from_module_context(cxt, 'counts')))
			
		labels = []
		files = []
		for fileinfo in self.resolve_input('bams', cxt):
			files.append(fileinfo.fullpath)
			labels.append(fileinfo.basename)
		fingerprint_args += [ '--bamfiles'] + files
		fingerprint_args += [ '--labels' ] + labels		

		cxt.log.write("Computing fingerprint for "+cxt.sample.name+" with DeepTools BamFingerprint......\n")
		cxt.log.write("DeepTools BamFingerprint version: "+self.tool_versions()['plotFingerprint']+"\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(fingerprint_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(fingerprint_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return out_files
	#end run()
#end class BamFingerprint