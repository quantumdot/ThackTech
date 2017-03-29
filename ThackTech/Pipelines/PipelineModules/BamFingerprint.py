import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class BamFingerprint(PipelineModule):
	
	def __init__(self, **kwargs):
		super(BamFingerprint, self).__init__('BamFingerprint', 'Bam Fingerprint', **kwargs)
		
		self._name_resolver('bams')
	#end __init__()
	
	def tool_versions(self):
		return {
			'plotFingerprint': subprocess.check_output("plotFingerprint --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'fingerprint')
		filetools.ensure_dir(dest_dir)
		fingerprint_args = [
			#'bamFingerprint',
			'plotFingerprint',	#renamed in deepTools 2.0
			'--numberOfProcessors', str(self.processors),
			'--plotFile',  os.path.join(dest_dir, cxt.sample.name+'.fingerprint.pdf'),
			'--plotFileFormat', 'pdf',
			'--outRawCounts', os.path.join(dest_dir, cxt.sample.name+'.fingerprint.rawcounts.txt')
		]
		labels = []
		files = []
		for fileinfo in self.resolve_input('bams', cxt):
			files.append(fileinfo.fullpath)
			labels.append(fileinfo.basename)
		fingerprint_args += [ '--bamfiles'] + files
		fingerprint_args += [ '--labels' ] + labels		

		cxt.log.write("Computing fingerprint for "+cxt.sample.name+" with DeepTools BamFingerprint......\n")
		cxt.log.write("DeepTools BamFingerprint version: "+subprocess.check_output(['plotFingerprint', '--version'], stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(fingerprint_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(fingerprint_args, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		return {
			'figure': os.path.join(dest_dir, cxt.sample.name+'.fingerprint.pdf'),
			'counts': os.path.join(dest_dir, cxt.sample.name+'.fingerprint.rawcounts.txt')
		}
	#end run()
#end class BamFingerprint