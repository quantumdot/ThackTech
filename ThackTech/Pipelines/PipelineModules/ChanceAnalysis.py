import os
import subprocess
import sys

from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ChanceAnalysis(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'Chance', 'IPStrength and spectrum analysis')
	#end __init__()

	def supported_types(self):
		return ['bed', 'bam', 'bampe', 'sam', 'mat', 'bowtie', 'tagalign']
	#end supported_types()
	
	def run(self, cxt):
		chance_path = '/home/josh/scripts/chance/run_chance_com.sh'
		dest_dir = os.path.join(cxt.sample.dest, 'chance')
		filetools.ensure_dir(dest_dir)
		output_files = {}
		
		shared_chance_args = [
			'-b', cxt.sample.genome.name,
			'-t', 'bam' if cxt.sample.format in ['bam', 'bampe'] else cxt.sample.format
		]

		ip_strength_args = [
			chance_path,
			'IPStrength',
			'--ipfile', cxt.sample.get_file('source', 'treatment'),
			'--ipcxt.sample', cxt.sample.name,
			'-o', os.path.join(dest_dir, cxt.sample.name+'.ipstrength.txt')
		] + shared_chance_args
		if cxt.sample.has_file('source', 'control'):
			ip_strength_args += [ '--inputfile', cxt.sample.get_file('source', 'control'), '--inputcxt.sample', cxt.sample.name+'_INPUT' ]

		cxt.log.write("Computing IP Strength with Chance......\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(ip_strength_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		self._run_subprocess(ip_strength_args, cwd=dest_dir, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		output_files['ip_strength'] = os.path.join(dest_dir, cxt.sample.name+'.ipstrength.txt')

		for label, file_path in cxt.sample.get_file_group('source').iteritems():
			chance_spectrum_args = [
				chance_path,
				'spectrum',
				'-o', os.path.join(dest_dir, cxt.sample.name+'.'+label+'.spectrum.txt'),
				'-s', cxt.sample.name+'_'+label,
				'-f', file_path
			] + shared_chance_args

			cxt.log.write("Computing cxt.sample spectrum for "+label+" with Chance......\n")
			cxt.log.write("\n..............................................\n")
			cxt.log.write(" ".join(chance_spectrum_args))
			cxt.log.write("\n..............................................\n")
			cxt.log.flush()
			self._run_subprocess(chance_spectrum_args, cwd=dest_dir, stderr=subprocess.STDOUT, stdout=cxt.log)
			
			output_files[label+'_spectrum'] = os.path.join(dest_dir, cxt.sample.name+'.'+label+'.spectrum.txt')
		
		return output_files
	#end run()
#end class ChanceAnalysis