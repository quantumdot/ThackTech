import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class BamFingerprint(PipelineModule):
	
	def __init__(self, **kwargs):
		super(BamFingerprint, self).__init__('BamFingerprint', 'Bam Fingerprint', **kwargs)
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
		if cxt.sample.has_file('source', 'control'):
			files.append('"%s"' % (cxt.sample.get_file('source', 'control'),))
			labels.append('Input')
		files.append('"%s"' % (cxt.sample.get_file('source', 'treatment'),))
		labels.append(cxt.sample.name)
		fingerprint_args += [ '--bamfiles', ' '.join(files) ]
		fingerprint_args += [ '--labels', ' '.join(labels) ]

		cxt.log.write("Computing fingerprint for "+cxt.sample.name+" with DeepTools BamFingerprint......\n")
		cxt.log.write("DeepTools BamFingerprint version: "+subprocess.check_output(['plotFingerprint', '--version'], stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(fingerprint_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(' '.join(fingerprint_args), shell=True, stderr=subprocess.STDOUT, stdout=cxt.log) #for some reason this doesnt seem to work when not run through a shell
																	#I suspect somehow the multiple args for --bamfiles are getting screwed up by python! eekk!

		
		return {
			'figure': os.path.join(dest_dir, cxt.sample.name+'.fingerprint.pdf'),
			'counts': os.path.join(dest_dir, cxt.sample.name+'.fingerprint.rawcounts.txt')
		}
	#end run()
#end class BamFingerprint