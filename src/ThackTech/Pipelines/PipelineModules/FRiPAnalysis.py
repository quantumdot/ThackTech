import os
import subprocess
import sys

import ThackTech.Common as Common
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule, ModuleParameter


class FRiPAnalysis(PipelineModule):
	
	def __init__(self):
		PipelineModule.__init__(self, 'FRiP', 'Fraction of Reads in Peaks')
		
		self._name_resolver('peaks')
	#end __init__()

	def supported_types(self):
		return ['bam', 'bampe']
	#end supported_types()
	
	def run(self, sample, logfile):
		dest_dir = os.path.join(sample.dest, 'FRiP')
		filetools.ensure_dir(dest_dir)
		with open(os.path.join(dest_dir, sample.name+'.FRiP.txt'), 'w') as results_file:
			#results_file.write("FRiP (Fraction of Reads in Peaks) for sample group %s\n" % (sample.name,))
			#results_file.write("================================================================================\n")
			results_file.write("Sample\tType\tMapped_Reads\tReads_In_Bed_Regions\tFRiP\n")
			
			for label, file_path in sample.get_file_group('source').iteritems():
				logfile.write("Computing FRiP for %s\n" % (label,))
				logfile.flush()
				peak_file = self.resolve_input('peaks', sample)
				if not os.path.isfile(peak_file):
					raise IOError('Unable to locate peaks file [%s]!' % (peak_file,))
				peak_type = os.path.splitext(peak_file)[1][1:]
				if peak_type == 'bed':
					cut_col = 6
				elif peak_type == 'broadPeak':
					cut_col = 10
				elif peak_type == 'narrowPeak':
					cut_col = 11
				else:
					logfile.write("Unable to find suitable intervals for %s\n" % (label,))
					return False

				total_mapped_reads   = float(subprocess.check_output('samtools flagstat "%s" | grep -Po "([0-9]+)(?= \+ [0-9]+ mapped)"' % (file_path,), shell=True))
				reads_in_bed_regions = float(subprocess.check_output('bedtools coverage -abam "%s" -b "%s" -counts | cut -f %d | paste -sd+ - | bc' % (file_path, peak_file, cut_col), shell=True))
				frip = reads_in_bed_regions / total_mapped_reads

				results_file.write("%s\t%s\t%d\t%d\t%f\n" % (sample.name, label, total_mapped_reads, reads_in_bed_regions, frip))
		return {
			'frip': os.path.join(dest_dir, sample.name+'.FRiP.txt')
		}
	#end run()
#end class FrIPAnalysis