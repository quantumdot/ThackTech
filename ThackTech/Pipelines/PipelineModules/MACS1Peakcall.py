import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class MACS1Peakcall(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MACS1', short_description='Peak Calling with MACS1')
		super_args.update(**kwargs)
		super(MACS1Peakcall, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		self.add_parameter(ModuleParameter('duplicates', str, 	'auto',	desc="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}."))
		self.add_parameter(ModuleParameter('bw', 		 int, 	300,	desc="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab."))
		self.add_parameter(ModuleParameter('sigout',	 str, 	'bdg',	desc="Output type for signal. Either 'wig' or 'bdg'."))
		self.add_parameter(ModuleParameter('tsize',		 int, 	None, nullable=True,	desc="Tag size. This will overide the auto detected tag size."))
		self.add_parameter(ModuleParameter('pvalue',	 float,	1e-5, desc="Pvalue cutoff for peak detection."))
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('treatments')
		self._name_resolver('controls')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'macs': self._call_output("macs --version 2>&1 | perl -ne 'if(m/macs\s+(.+)$/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):

		macs_args = [
			'macs',
			'--name', cxt.sample.name,
			'--gsize', str(cxt.sample.genome.gsize),
			#'--format', self.get_file_format(treatment),
			'--keep-dup', self.get_parameter_value_as_string('duplicates'),
			'--bw', self.get_parameter_value_as_string('bw'),
			'--verbose', '2', 
			'--single-profile', 
			'--diag',
			('--bdg' if self.get_parameter_value_as_string('sigout') == 'bdg' else '--wig'),
		]
		
		treatments = self.resolve_input('treatments', cxt)
		macs_args.extend(['--treatment'] + [f.fullpath for f in treatments])
		
		controls = self.resolve_input('controls', cxt)
		if controls is not None and len(controls) > 0:
			macs_args.extend(['--control'] + [f.fullpath for f in controls])

			
		cxt.log.write("Performing peak calling with MACS......\n")
		cxt.log.write("-> "+subprocess.check_output(['macs', '--version'], stderr=subprocess.STDOUT)+"")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(macs_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()	
		self._run_subprocess(macs_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		if self.get_parameter_value_as_string('sigout') == 'bdg':
			signal_output_ext = '.bdg.gz'
			signal_folder = '_MACS_bedGraph'
		else:
			signal_output_ext = '.wig.gz'
			signal_folder = '_MACS_wiggle'
		output_files = {
			'peaks_xls':		os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.xls'),
			'neg_peaks_xls':	os.path.join(cxt.sample.dest, cxt.sample.name+'_negative_peaks.xls'),
			'peaks':			os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.bed'),
			'summits':			os.path.join(cxt.sample.dest, cxt.sample.name+'_summits.bed'),
			'treatment_signal':	os.path.join(cxt.sample.dest, cxt.sample.name+signal_folder, 'treat', cxt.sample.name+'_treat_afterfiting_all'+signal_output_ext)
		}
		
		if os.path.isfile(os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')):
			output_files['model_Rscript'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')
			cxt.log.write('Generating MACS model figure.....\n')
			cxt.log.flush()
			output_files['model_figure'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.pdf')
			self._run_subprocess(['Rscript', '--vanilla', output_files['model_Rscript']], stderr=subprocess.STDOUT, stdout=cxt.log, cwd=cxt.sample.dest)
		else:
			cxt.log.write('Skipping MACS model figure generation because model R script is not present.....\n')
			cxt.log.flush()

		if controls is not None and len(controls) > 0:
			output_files['control_signal'] = os.path.join(cxt.sample.dest, cxt.sample.name+signal_folder, 'control', cxt.sample.name+'_control_afterfiting_all'+signal_output_ext)
		
		return output_files
	#end run()
	
	def get_file_format(self, fileinfo):
		ext = fileinfo.ext.lower()
		if ext == '.bam':
			return 'BAM'
		elif ext == '.sam':
			return 'SAM'
		elif ext == '.bed':
			return 'BED'
		elif ext == '.bowtie':
			return 'BOWTIE'
		else:
			return 'AUTO'
	#end get_file_format()
#end class MACS1Peakcall
