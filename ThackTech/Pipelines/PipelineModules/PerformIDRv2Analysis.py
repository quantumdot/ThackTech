import itertools
import math
import os
import subprocess
import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class PerformIDRv2Analysis(PipelineModule):
	
	def __init__(self, **kwargs):
		super(PerformIDRv2Analysis, self).__init__('IDRv2', 'Perform IDRv2 Analysis', **kwargs)
		self._name_resolver('primary_replicates')
		self._name_resolver('pseudo_replicates')
		self._name_resolver('pooled_pseudo_replicates')
		
		#self.add_parameter('peak_half_width', '-1')
		#self.add_parameter('min_overlap_ratio', '0')
		self.add_parameter(ModuleParameter('ranking_measure', str, 'p.value'))
		
		self.add_parameter(ModuleParameter('primary_replicates_IDR_threshold', float, 0.01))
		self.add_parameter(ModuleParameter('pseudo_replicates_IDR_threshold', float, 0.01))
		self.add_parameter(ModuleParameter('pooled_pseudo_replicates_IDR_threshold', float, 0.0025))
		
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		curr_dest = cxt.sample.dest
		cxt.sample.dest = os.path.join(cxt.sample.dest, 'idr')
		filetools.ensure_dir(cxt.sample.dest)
		output_files = {}
		
		consistancy_output = os.path.join(cxt.sample.dest, cxt.sample.name+'_idr_npeaks.txt')
		with open(consistancy_output, 'a') as count_file:
			count_file.write('Group\tComparison\tThreshold\tNumPeaks\n')
			
			#IDR for primary replicates
			output = self._run_batch_consistency_analysis('primary_replicates', self.get_parameter_value('primary_replicates_IDR_threshold'), cxt)
			output_files.update(output)
			for key in output:
				if key.endswith('_N'):
					count_file.write('%s\t%s\t%f\t%d\n' % ('primary_replicates', key.replace('_N', ''), self.get_parameter_value('primary_replicates_IDR_threshold'), output[key]))
			
			#IDR for pseudo self-replicates
			output = self._run_batch_consistency_analysis('pseudo_replicates',  self.get_parameter_value('pseudo_replicates_IDR_threshold'), cxt)
			output_files.update(output)
			for key in output:
				if key.endswith('_N'):
					count_file.write('%s\t%s\t%f\t%d\n' % ('pseudo_replicates', key.replace('_N', ''), self.get_parameter_value('pseudo_replicates_IDR_threshold'), output[key]))
			
			#IDR for pseudo pooled-replicates
			output = self._run_batch_consistency_analysis('pooled_pseudo_replicates', self.get_parameter_value('pooled_pseudo_replicates_IDR_threshold'), cxt)
			output_files.update(output)
			for key in output:
				if key.endswith('_N'):
					count_file.write('%s\t%s\t%f\t%d\n' % ('pooled_pseudo_replicates', key.replace('_N', ''), self.get_parameter_value('pooled_pseudo_replicates_IDR_threshold'), output[key]))
		
		output_files['consistancy_output'] = consistancy_output
		cxt.sample.dest = curr_dest #replace the cxt.sample destination
		return output_files
	#end run()
	
	def _run_batch_consistency_analysis(self, replicate_type, idr_threshold, cxt):
		output_files = {}
		replicate_combinations = list(itertools.combinations(self.resolve_input(replicate_type, cxt.sample), 2))
		output_prefixes = []
		for pair in replicate_combinations:
			rep1_bn = filetools.basename_noext(pair[0])
			rep2_bn = filetools.basename_noext(pair[1])
			
			peak_type = os.path.splitext(pair[0])[1][1:]
			output_name = rep1_bn+'_VS_'+rep2_bn
			pooled_common_peaks_IDR_filename = os.path.join(cxt.sample.dest, output_name+'.pooled_common_IDRv2.'+peak_type)
			idr_log_output = os.path.join(cxt.sample.dest, output_name+'.IDR.log')
			idr_cmd = [ 
				'idr',
				'--verbose',
				'--plot',
				'--input-file-type', peak_type,
				'--rank', self.get_parameter_value('ranking_measure'),
				'--output-file', pooled_common_peaks_IDR_filename,
				'--log-output-file', idr_log_output,
				'--soft-idr-threshold', str(idr_threshold),
				'--cxt.samples', pair[0], pair[1]
			]
			

			cxt.log.write('-> Performing IDR analysis on %s VS %s\n' % (rep1_bn, rep2_bn))
			cxt.log.write("-> "+subprocess.check_output(['idr', '--version'], stderr=subprocess.STDOUT))
			cxt.log.write("..............................................\n")
			cxt.log.write(" ".join(idr_cmd))
			cxt.log.write("\n..............................................\n")
			cxt.log.flush()
			proc = subprocess.Popen(idr_cmd, stderr=subprocess.STDOUT, stdout=cxt.log)
			proc.communicate()

			output_files[output_name + '_pooled_common_peaks_IDR'] = pooled_common_peaks_IDR_filename
			output_files[output_name + '_EM_parameters_log'] = idr_log_output
			output_files[output_name + '_IDR_plot']	= os.path.join(cxt.sample.dest, output_name+'.png')
			
		
		
			awk_string = r"""awk 'BEGIN{OFS="\t"} $12>=%2.2f {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'""" % (-math.log10(idr_threshold))
			final_IDR_thresholded_filename =  os.path.join(cxt.sample.dest, output_name + '.IDR%2.2f.%s' % (idr_threshold, peak_type))
			Common.run_pipe([
				'cat %s' % (pooled_common_peaks_IDR_filename),
				awk_string,
				'sort -k7n,7n'
				#'gzip -c'
			], final_IDR_thresholded_filename)

			npeaks_pass_filename = os.path.join(cxt.sample.dest, output_name + '-npeaks-aboveIDR.txt')
			wc_output = subprocess.check_output(['wc', '-l', final_IDR_thresholded_filename])
			with open(npeaks_pass_filename, 'w') as fh:
				fh.write(wc_output)
			line_count = wc_output.split()[0]
			n_peaks = int(line_count)
			
			output_files[output_name + '_final_IDR_thresholded_peaks'] = final_IDR_thresholded_filename
			output_files[output_name + '_npeaks_pass_IDR'] = npeaks_pass_filename
			output_files[output_name + '_N'] = n_peaks

		
		return output_files
	#end _run_batch_consistency_analysis()
#end class PerformIDRAnalysis