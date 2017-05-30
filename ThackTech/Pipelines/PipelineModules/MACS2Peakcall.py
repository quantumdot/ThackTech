import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MACS2Peakcall(PipelineModule):
	
	def __init__(self, **kwargs):
		super(MACS2Peakcall, self).__init__('MACS2', 'Peak calling with MACS2', **kwargs)
		
		self.add_parameter(ModuleParameter('duplicates', str, 	'auto',	desc="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}."))
		self.add_parameter(ModuleParameter('bw', 		 int, 	300,	desc="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab."))
		self.add_parameter(ModuleParameter('sigout', 	str, 	'bdg',	desc="Output type for signal. Either 'wig' or 'bdg'."))
		
		self._name_resolver('treatments')
		self._name_resolver('controls')
	#end __init__()

	def supported_types(self):
		return ["bed", "eland", "elandmulti", "elandmultipet", "elandexport", "sam", "bam", "bampe", "bowtie"]
	#end supported_types()
	
	def tool_versions(self):
		return {
			'macs2': self._call_output("macs2 --version 2>&1 | perl -ne 'if(m/macs2\s+(.+)$/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		macs_args = [
			'macs2',
			'callpeak',
			'--verbose', '3',
			'--name', cxt.sample.name,
			'--gsize', str(cxt.sample.genome.gsize),
			#'--format', cxt.sample.format.upper(),
			'--format', 'BAM',
			'--keep-dup', self.get_parameter_value_as_string('duplicates'),
			'--bw', self.get_parameter_value_as_string('bw'),
			'--cutoff-analysis',
			'--bdg'
		]

		treatments = self.resolve_input('treatments', cxt)
		macs_args.extend(['--treatment'] + [f.fullpath for f in treatments])
		
		controls = self.resolve_input('controls', cxt)
		if len(controls) > 0:
			macs_args.extend(['--control'] + [f.fullpath for f in controls])

		
			
		if cxt.sample.has_attribute('broad') and cxt.sample.get_attribute('broad'):
			macs_args.append('--broad')
			macs_args += [ '--broad-cutoff', '0.1' ]
			
		cxt.log.write("Performing peak calling with MACS......\n")
		#cxt.log.write("-> "+subprocess.check_output('/bin/bash -c "source /home/josh/scripts/macs2dev/bin/activate && macs2 --version"', shell=True, stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("-> "+self._call_output('macs2 --version', shell=True, stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(macs_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()		

		self._run_subprocess(macs_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		output_files = {
			'cutoff_analysis':	os.path.join(cxt.sample.dest, cxt.sample.name+'_cutoff_analysis.txt'),
			'treatment_signal':	os.path.join(cxt.sample.dest, cxt.sample.name+'_treat_pileup.bdg'),
			'peaks_xls':		os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.xls'),
		}
		if os.path.isfile(os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')):
			output_files['model_Rscript'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')
			
		if len(controls) > 0:
			output_files['control_signal'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_control_lambda.bdg')
			
		if cxt.sample.has_attribute('broad') and cxt.sample.get_attribute('broad'):
			output_files['broad_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.broadPeak')
			output_files['gapped_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.gappedPeak')
		else:
			output_files['narrow_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.narrowPeak')
			output_files['summits'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_summits.bed')
			
		return output_files
	#end run()
#end class MACS2Peakcall