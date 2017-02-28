import os
import subprocess
import sys

from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MACS2Peakcall(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'MACS2', 'Peak calling with MACS2')
		
		self.add_parameter(ModuleParameter('duplicates', str, 	'auto',	desc="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}."))
		self.add_parameter(ModuleParameter('bw', 		 int, 	300,	desc="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab."))
		self.add_parameter(ModuleParameter('sigout', 	str, 	'bdg',	desc="Output type for signal. Either 'wig' or 'bdg'."))
	#end __init__()

	def supported_types(self):
		return ["bed", "eland", "elandmulti", "elandmultipet", "elandexport", "sam", "bam", "bampe", "bowtie"]
	#end supported_types()
	
	def tool_versions(self):
		return {
			'macs2': subprocess.check_output("macs2 --version 2>&1 | perl -ne 'if(m/macs2\s+(.+)$/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, sample, logfile):
		macs_args = [
			'macs2',
			'callpeak',
			'--verbose', '3',
			'--name', sample.name,
			'--gsize', sample.genome.gsize,
			'--format', sample.format.upper(),
			'--keep-dup', self.get_parameter_value_as_string('duplicates'),
			'--bw', self.get_parameter_value_as_string('bw'),
			'--cutoff-analysis',
			'--bdg',
			'--treatment', sample.get_file('source', 'treatment')
		]

		if sample.has_file('source', 'control'):
			macs_args += [ '--control', sample.get_file('source', 'control') ]
			
		if sample.has_attribute('broad') and sample.get_attribute('broad'):
			macs_args.append('--broad')
			macs_args += [ '--broad-cutoff', '0.1' ]
			
		logfile.write("Performing peak calling with MACS......\n")
		#logfile.write("-> "+subprocess.check_output('/bin/bash -c "source /home/josh/scripts/macs2dev/bin/activate && macs2 --version"', shell=True, stderr=subprocess.STDOUT)+"\n")
		logfile.write("-> "+subprocess.check_output('macs2 --version', shell=True, stderr=subprocess.STDOUT)+"\n")
		logfile.write("\n..............................................\n")
		logfile.write(" ".join(macs_args))
		logfile.write("\n..............................................\n")
		logfile.flush()		

		self._run_subprocess(macs_args, cwd=sample.dest, stderr=subprocess.STDOUT, stdout=logfile)
		
		output_files = {
			'cutoff_analysis':	os.path.join(sample.dest, sample.name+'_cutoff_analysis.txt'),
			'treatment_signal':	os.path.join(sample.dest, sample.name+'_treat_pileup.bdg'),
			'peaks_xls':		os.path.join(sample.dest, sample.name+'_peaks.xls'),
		}
		if os.path.isfile(os.path.join(sample.dest, sample.name+'_model.r')):
			output_files['model_Rscript'] = os.path.join(sample.dest, sample.name+'_model.r')
		if sample.has_file('source', 'control'):
			output_files['control_signal'] = os.path.join(sample.dest, sample.name+'_control_lambda.bdg')
		if sample.has_attribute('broad') and sample.get_attribute('broad'):
			output_files['broad_peaks'] = os.path.join(sample.dest, sample.name+'_peaks.broadPeak')
			output_files['gapped_peaks'] = os.path.join(sample.dest, sample.name+'_peaks.gappedPeak')
		else:
			output_files['narrow_peaks'] = os.path.join(sample.dest, sample.name+'_peaks.narrowPeak')
			output_files['summits'] = os.path.join(sample.dest, sample.name+'_summits.bed')
			
		return output_files
	#end run()
#end class MACS2Peakcall