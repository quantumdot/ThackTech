import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class ChanceAnalysis(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'Chance', 'IPStrength and spectrum analysis')
	#end __init__()

	def supported_types(self):
		return ['bed', 'bam', 'bampe', 'sam', 'mat', 'bowtie', 'tagalign']
	#end supported_types()
	
	def run(self, sample, logfile):
		chance_path = '/home/josh/scripts/chance/run_chance_com.sh'
		dest_dir = os.path.join(sample.dest, 'chance')
		Common.ensure_dir(dest_dir)
		output_files = {}
		
		shared_chance_args = [
			'-b', sample.genome.name,
			'-t', 'bam' if sample.format in ['bam', 'bampe'] else sample.format
		]

		ip_strength_args = [
			chance_path,
			'IPStrength',
			'--ipfile', sample.get_file('source', 'treatment'),
			'--ipsample', sample.name,
			'-o', os.path.join(dest_dir, sample.name+'.ipstrength.txt')
		] + shared_chance_args
		if sample.has_file('source', 'control'):
			ip_strength_args += [ '--inputfile', sample.get_file('source', 'control'), '--inputsample', sample.name+'_INPUT' ]

		logfile.write("Computing IP Strength with Chance......\n")
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(ip_strength_args))
		logfile.write("\n..............................................\n")
		logfile.flush()
		self._run_subprocess(ip_strength_args, cwd=dest_dir, stderr=subprocess.STDOUT, stdout=logfile)
		
		output_files['ip_strength'] = os.path.join(dest_dir, sample.name+'.ipstrength.txt')

		for label, file_path in sample.get_file_group('source').iteritems():
			chance_spectrum_args = [
				chance_path,
				'spectrum',
				'-o', os.path.join(dest_dir, sample.name+'.'+label+'.spectrum.txt'),
				'-s', sample.name+'_'+label,
				'-f', file_path
			] + shared_chance_args

			logfile.write("Computing sample spectrum for "+label+" with Chance......\n")
			logfile.write("\n..............................................\n")
			logfile.write(" ".join(chance_spectrum_args))
			logfile.write("\n..............................................\n")
			logfile.flush()
			self._run_subprocess(chance_spectrum_args, cwd=dest_dir, stderr=subprocess.STDOUT, stdout=logfile)
			
			output_files[label+'_spectrum'] = os.path.join(dest_dir, sample.name+'.'+label+'.spectrum.txt')
		
		return output_files
	#end run()
#end class ChanceAnalysis