import os
import subprocess
from ThackTech.Pipelines import PipelineModule, ModuleParameter
from ThackTech import aligntools


class MACS2Peakcall(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='MACS2', short_description='Peak calling with MACS2')
		super_args.update(**kwargs)
		super(MACS2Peakcall, self).__init__(**super_args)
	#end __init__()
	
	def __declare_parameters(self):
		#input arguments
		self.add_parameter(ModuleParameter('duplicates', str, 'auto', desc="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}."))
		
		#output arguments
		self.add_parameter(ModuleParameter('save_bedgraph', bool, True, desc="Whether to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file."))
		self.add_parameter(ModuleParameter('bdg_trackline', bool, False, desc="Include trackline with bedGraph files."))
		self.add_parameter(ModuleParameter('spmr', bool, True, desc="Save signal as signal per mission reads."))
		self.add_parameter(ModuleParameter('verbosity', int, 3, desc="Set verbose level of runtime message. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show debug messages."))
	
		#shift model arguments
		self.add_parameter(ModuleParameter('tag_size', int, None, nullable=True, desc="Size of sequencing tags. If None, then MACS will determine automatically."))
		self.add_parameter(ModuleParameter('bandwith', int, 300, desc="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab."))
		#--mfold
		#--fix-bimodal
		#--nomodel
		#--shift
		#--extsize
		
		#peak calling arguments
		self.add_parameter(ModuleParameter('pvalue', float, None, nullable=True, desc="Pvalue cutoff for peak detection."))
		self.add_parameter(ModuleParameter('qvalue', float, 0.05, desc="Minimum FDR (q-value) cutoff for peak detection."))
		self.add_parameter(ModuleParameter('broad_cutoff', float, 0.1, desc="Cutoff for broad region."))
		#--to-large
		#--ratio
		#--down-sample
		#--seed
		#--tempdir
		#--nolambda
		#--slocal
		#--llocal
		self.add_parameter(ModuleParameter('cutoff_analysis', bool, True, desc="Perform cutoff analysis."))
		
		#post-processing option
		#--call-summits
		#--fe-cutoff
		
		#end __declare_parameters()
	
	def __declare_resolvers(self):
		self._name_resolver('treatments')
		self._name_resolver('controls')
	#end __declare_resolvers()
	
	def tool_versions(self):
		return {
			'macs2': self._call_output("macs2 --version 2>&1 | perl -ne 'if(m/macs2\s+(.+)$/){ print $1 }'", shell=True, stderr=subprocess.STDOUT)
		}
	#end tool_versions()
	
	def run(self, cxt):
		output_files = {}
		
		macs_args = [
			'macs2',
			'callpeak',
			'--verbose', self.get_parameter_value_as_string('verbosity'),
			'--name', cxt.sample.name,
			'--gsize', str(cxt.sample.genome.gsize),
			'--outdir', cxt.sample.dest,
			'--keep-dup', self.get_parameter_value_as_string('duplicates'),
			'--bw', self.get_parameter_value_as_string('bw'),
		]
		
		if self.get_parameter_value('save_bedgraph'):
			macs_args.append('--bdg')
			if self.get_parameter_value('spmr'):
				macs_args.append('--SPMR')
			if self.get_parameter_value('bdg_trackline'):
				macs_args.append('--trackline')
		
		if self.get_parameter_value('cutoff_analysis'):
			macs_args.append('--cutoff-analysis')
			output_files['cutoff_analysis'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_cutoff_analysis.txt')
		
		if self.get_parameter_value('tag_size') is not None:
			macs_args.extend(['--tsize', self.get_parameter_value_as_string('tag_size')])
			
		if self.get_parameter_value('pvalue') is not None:
			macs_args.extend(['--pvalue', self.get_parameter_value_as_string('pvalue')])
		else:
			macs_args.extend(['--qvalue', self.get_parameter_value_as_string('qvalue')])
			
		if cxt.sample.has_attribute('broad') and cxt.sample.get_attribute('broad'):
			macs_args.extend(['--broad', '--broad-cutoff', '0.1' ])

		formats = []
		treatments = self.resolve_input('treatments', cxt)
		macs_args.extend(['--treatment'] + [f.fullpath for f in treatments])
		formats.extend([self.get_file_format(f) for f in treatments])
		
		controls = self.resolve_input('controls', cxt)
		if controls is not None and len(controls) > 0:
			macs_args.extend(['--control'] + [f.fullpath for f in controls])
			formats.extend([self.get_file_format(f) for f in treatments])

		formats = set(formats)
		if len(formats) > 1:
			#macs reccomends using 'auto' if multiple formats are given
			macs_args.extend(['--format', 'AUTO'])
		else:
			macs_args.extend(['--format', list(formats)[0]])
		
			
		
			
		cxt.log.write("Performing peak calling with MACS......\n")
		#cxt.log.write("-> "+subprocess.check_output('/bin/bash -c "source /home/josh/scripts/macs2dev/bin/activate && macs2 --version"', shell=True, stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("-> "+self._call_output('macs2 --version', shell=True, stderr=subprocess.STDOUT)+"\n")
		cxt.log.write("\n..............................................\n")
		cxt.log.write(" ".join(macs_args))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()		

		self._run_subprocess(macs_args, cwd=cxt.sample.dest, stderr=subprocess.STDOUT, stdout=cxt.log)
		
		
		
		if os.path.isfile(os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')):
			output_files['model_Rscript'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.r')
			cxt.log.write('Generating MACS model figure.....\n')
			cxt.log.flush()
			output_files['model_figure'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_model.pdf')
			self._run_subprocess(['Rscript', '--vanilla', output_files['model_Rscript']], stderr=subprocess.STDOUT, stdout=cxt.log, cwd=cxt.sample.dest)
		else:
			cxt.log.write('Skipping MACS model figure generation because model R script is not present.....\n')
			cxt.log.flush()
		
		
		
		output_files['treatment_signal'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_treat_pileup.bdg')
		output_files['peaks_xls'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.xls')
		
		if controls is not None and len(controls) > 0:
			output_files['control_signal'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_control_lambda.bdg')
			
		if cxt.sample.has_attribute('broad') and cxt.sample.get_attribute('broad'):
			output_files['broad_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.broadPeak')
			output_files['gapped_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.gappedPeak')
		else:
			output_files['narrow_peaks'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_peaks.narrowPeak')
			output_files['summits'] = os.path.join(cxt.sample.dest, cxt.sample.name+'_summits.bed')
			
		return output_files
	#end run()
	
	
	def get_file_format(self, fileinfo):
		ext = fileinfo.ext.lower()
		if ext == '.bam':
			if aligntools.is_bam_PE(fileinfo.fullpath):
				return 'BAMPE'
			else:
				return 'BAM'
		elif ext == '.sam':
			return 'SAM'
		elif ext == '.bed':
			return 'BED'
		elif ext == '.bedpe':
			return 'BEDPE'
		elif ext == '.bowtie':
			return 'BOWTIE'
		else:
			return 'AUTO'
	#end get_file_format()
#end class MACS2Peakcall