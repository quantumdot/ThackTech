import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MACS1Peakcall(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'MACS1', 'Peak Calling with MACS1')
		
		self.add_parameter(ModuleParameter('duplicates', str, 	'auto',	desc="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}."))
		self.add_parameter(ModuleParameter('bw', 		 int, 	300,	desc="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab."))
		self.add_parameter(ModuleParameter('sigout', 	str, 	'bdg',	desc="Output type for signal. Either 'wig' or 'bdg'."))
	#end __init__()
	
	def tool_versions(self):
		return {
			'macs': subprocess.check_output("macs --version 2>&1 | perl -ne 'if(m/macs\s+(.+)$/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, sample, logfile):
		format = sample.format.upper()
		if format == "BAMPE":
			format = "BAM"

		macs_args = [
			'macs',
			'--name', sample.name,
			'--gsize', sample.genome.gsize,
			'--format', format,
			'--keep-dup', self.get_parameter_value_as_string('duplicates'),
			'--bw', self.get_parameter_value_as_string('bw'),
			'--verbose', '2', 
			'--single-profile', 
			'--diag',
			('--bdg' if self.get_parameter_value_as_string('sigout') == 'bdg' else '--wig'),
			'--treatment', sample.get_file('source', 'treatment')
		]

		if sample.has_file('source', 'control'):
			macs_args += [ '--control', sample.get_file('source', 'control') ]
			
		logfile.write("Performing peak calling with MACS......\n")
		logfile.write("-> "+subprocess.check_output(['macs', '--version'], stderr=subprocess.STDOUT)+"")
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(macs_args))
		logfile.write("\n..............................................\n")
		logfile.flush()	
		self._run_subprocess(macs_args, cwd=sample.dest, stderr=subprocess.STDOUT, stdout=logfile)
		
		if self.get_parameter_value_as_string('sigout') == 'bdg':
			signal_output_ext = '.bdg.gz'
			signal_folder = '_MACS_bedGraph'
		else:
			signal_output_ext = '.wig.gz'
			signal_folder = '_MACS_wiggle'
		output_files = {
			'peaks_xls':		os.path.join(sample.dest, sample.name+'_peaks.xls'),
			'neg_peaks_xls':	os.path.join(sample.dest, sample.name+'_negative_peaks.xls'),
			'peaks':			os.path.join(sample.dest, sample.name+'_peaks.bed'),
			'summits':			os.path.join(sample.dest, sample.name+'_summits.bed'),
			'treatment_signal':	os.path.join(sample.dest, sample.name+signal_folder, 'treat', sample.name+'_treat_afterfiting_all'+signal_output_ext)
		}
		if os.path.isfile(os.path.join(sample.dest, sample.name+'_model.r')):
			output_files['model_Rscript'] = os.path.join(sample.dest, sample.name+'_model.r')
		if sample.has_file('source', 'control'):
			output_files['control_signal'] = os.path.join(sample.dest, sample.name+signal_folder, 'control', sample.name+'_control_afterfiting_all'+signal_output_ext)
		return output_files
	#end run()
#end class MACS1Peakcall
