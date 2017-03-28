import itertools
import os
import random
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter
import pysam


class PerformIDRAnalysis(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'IDR', 'Perform IDR Analysis')
		self._name_resolver('primary_replicates')
		self._name_resolver('pseudo_replicates')
		self._name_resolver('pooled_pseudo_replicates')
		
		self.add_parameter(ModuleParameter('peak_half_width', int, '-1'))
		self.add_parameter(ModuleParameter('min_overlap_ratio', float, '0'))
		self.add_parameter(ModuleParameter('ranking_measure', str, 'p.value'))
		
		self.add_parameter(ModuleParameter('primary_replicates_IDR_threshold', float, 0.01))
		self.add_parameter(ModuleParameter('pseudo_replicates_IDR_threshold', float, 0.01))
		self.add_parameter(ModuleParameter('pooled_pseudo_replicates_IDR_threshold', float, 0.0025))
		
	#end __init__()

	def supported_types(self):
		return None
	#end supported_types()
	
	def run(self, cxt):
		dest_dir = os.path.join(cxt.sample.dest, 'idr')
		filetools.ensure_dir(dest_dir)
		output_files = {}
		
		consistancy_output = os.path.join(dest_dir, cxt.sample.name+'_idr_npeaks.txt')
		with open(consistancy_output, 'a') as count_file:
			count_file.write('Group\tComparison\tThreshold\tNumPeaks\n')
			output = self._run_batch_consistency_analysis('primary_replicates', cxt.sample, dest_dir, cxt.sample.has_attribute('broad'), cxt.log)
			output_files.update(output)
			for key in output:
				if key.endswith('_overlapping_peaks'):
					n = self._get_num_consistant_peaks(output_files[key], self.get_parameter_value('primary_replicates_IDR_threshold'))
					count_file.write('%s\t%f\t%d\n' % ('primary_replicates', key.replace('_overlapping_peaks', ''), self.get_parameter_value('primary_replicates_IDR_threshold'), n))
			
			output = self._run_batch_consistency_analysis('pseudo_replicates',  cxt.sample, dest_dir, cxt.sample.has_attribute('broad'), cxt.log)
			output_files.update(output)
			for key in output:
				if key.endswith('_overlapping_peaks'):
					n = self._get_num_consistant_peaks(output_files[key], self.get_parameter_value('pseudo_replicates_IDR_threshold'))
					count_file.write('%s\t%f\t%d\n' % ('pseudo_replicates', key.replace('_overlapping_peaks', ''), self.get_parameter_value('pseudo_replicates_IDR_threshold'), n))
					
			output = self._run_batch_consistency_analysis('pooled_pseudo_replicates',  cxt.sample, dest_dir, cxt.sample.has_attribute('broad'), cxt.log)
			output_files.update(output)
			for key in output:
				if key.endswith('_overlapping_peaks'):
					n = self._get_num_consistant_peaks(output_files[key], self.get_parameter_value('pooled_pseudo_replicates_IDR_threshold'))
					count_file.write('%s\t%f\t%d\n' % ('pooled_pseudo_replicates', key.replace('_overlapping_peaks', ''), self.get_parameter_value('pooled_pseudo_replicates_IDR_threshold'), n))
		
		output_files['consistancy_output'] = consistancy_output
		return output_files
	#end run()
	
	def _get_num_consistant_peaks(self, overlap_peak_file, threshold):
		result = subprocess.check_output("awk '$11 <= "+str(threshold)+" {print $0}' "+overlap_peak_file+" | wc -l", shell=True, stderr=subprocess.STDOUT)
		return int(result)
	#end _get_num_consistant_peaks()
	
	def _run_batch_consistency_analysis(self, replicate_type, cxt.sample, dest, broad, cxt.log):
		output_files = {}
		replicate_combinations = list(itertools.combinations(self.resolve_input(replicate_type, cxt.sample), 2))
		output_prefixes = []
		for pair in replicate_combinations:
			rep1_bn = Common.basename_noext(pair[0])
			rep2_bn = Common.basename_noext(pair[1])
			output_name = rep1_bn+'_VS_'+rep2_bn
			outprefix = os.path.join(dest, output_name)
			output_prefixes.append(outprefix)
			idr_cmd = [ 
				'Rscript', 'batch-consistency-analysis.r',
				pair[0], pair[1],
				self.get_parameter_value('peak_half_width'),
				outprefix,
				self.get_parameter_value('peak_half_width'),				
				'T' if broad else 'F',
				self.get_parameter_value('ranking_measure')
			]
			cxt.log.write('-> Performing IDR analysis on %s VS %s\n' % (rep1_bn, rep2_bn))
			cxt.log.write("..............................................\n")
			cxt.log.write(" ".join(idr_cmd))
			cxt.log.write("\n..............................................\n")
			cxt.log.flush()
			proc = subprocess.Popen(idr_cmd, cwd='/home/josh/scripts/idrCode/', stderr=subprocess.STDOUT, stdout=cxt.log)
			proc.communicate()

			output_files[output_name + '_EM_fitting_output'] 	= outprefix + '-em.sav'
			output_files[output_name + '_emperical_curve_data'] = outprefix + '-uri.sav'
			output_files[output_name + '_params_logs'] 			= outprefix + '-Rout.txt'
			output_files[output_name + '_npeaks_aboveIDR']		= outprefix + '-npeaks-aboveIDR.txt'
			output_files[output_name + '_overlapping_peaks']	= outprefix + '-overlapped-peaks.txt'
		
		###
		# Generate IDR plots
		###
		idr_plot_cmd = [
			'Rscript', 'batch-consistency-plot.r',
			str(len(replicate_combinations)),
			os.path.join(dest, replicate_type)
		] + output_prefixes
		cxt.log.write('-> Generating IDR plots for %s\n' % (replicate_type,))
		cxt.log.write("..............................................\n")
		cxt.log.write(" ".join(idr_plot_cmd))
		cxt.log.write("\n..............................................\n")
		cxt.log.flush()
		proc = subprocess.Popen(idr_cmd, cwd='/home/josh/scripts/idrCode/', stderr=subprocess.STDOUT, stdout=cxt.log)
		proc.communicate()
		output_files[replicate_type + '_IDR_plot'] = os.path.join(dest, replicate_type+'-plot.ps')
		
		return output_files
	#end _run_batch_consistency_analysis()
#end class PerformIDRAnalysis