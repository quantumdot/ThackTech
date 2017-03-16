import os
import subprocess
import sys

from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class BamFingerprint(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'bamFingerprint', 'Bam Fingerprint')
	#end __init__()
	
	def tool_versions(self):
		return {
			'plotFingerprint': subprocess.check_output("plotFingerprint --version 2>&1 | perl -ne 'if(m/([\d\.]+)/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, sample, logfile):
		dest_dir = os.path.join(sample.dest, 'fingerprint')
		filetools.ensure_dir(dest_dir)
		fingerprint_args = [
			#'bamFingerprint',
			'plotFingerprint',	#renamed in deepTools 2.0
			'--numberOfProcessors', str(self.processors),
			'--plotFile',  os.path.join(dest_dir, sample.name+'.fingerprint.pdf'),
			'--plotFileFormat', 'pdf',
			'--outRawCounts', os.path.join(dest_dir, sample.name+'.fingerprint.rawcounts.txt')
		]
		labels = []
		files = []
		if sample.has_file('source', 'control'):
			files.append('"%s"' % (sample.get_file('source', 'control'),))
			labels.append('Input')
		files.append('"%s"' % (sample.get_file('source', 'treatment'),))
		labels.append(sample.name)
		fingerprint_args += [ '--bamfiles', ' '.join(files) ]
		fingerprint_args += [ '--labels', ' '.join(labels) ]

		logfile.write("Computing fingerprint for "+sample.name+" with DeepTools BamFingerprint......\n")
		logfile.write("DeepTools BamFingerprint version: "+subprocess.check_output(['plotFingerprint', '--version'], stderr=subprocess.STDOUT)+"\n")
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(fingerprint_args))
		logfile.write("\n..............................................\n")
		logfile.flush()
		self._run_subprocess(' '.join(fingerprint_args), shell=True, stderr=subprocess.STDOUT, stdout=logfile) #for some reason this doesnt seem to work when not run through a shell
																	#I suspect somehow the multiple args for --bamfiles are getting screwed up by python! eekk!

		
		return {
			'figure': os.path.join(dest_dir, sample.name+'.fingerprint.pdf'),
			'counts': os.path.join(dest_dir, sample.name+'.fingerprint.rawcounts.txt')
		}
	#end run()
#end class BamFingerprint