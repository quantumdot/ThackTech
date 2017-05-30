import os
import subprocess
from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule


class FRiPAnalysis(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='FRiP', short_description='Fraction of Reads in Peaks')
		super_args.update(**kwargs)
		super(FRiPAnalysis, self).__init__(**super_args)
		
		self._name_resolver('bed')
		self._name_resolver('bams')
	#end __init__()
	
	def run(self, cxt):
		
		peak_file = self.resolve_input('bed', cxt)
		if not peak_file.isfile:
			raise IOError('Unable to locate peak file {}!'.format(peak_file.fullpath))
		peak_type = peak_file.ext.lower()
		if peak_type == '.bed':
			cut_col = 6
		elif peak_type == '.broadpeak':
			cut_col = 10
		elif peak_type == '.narrowpeak':
			cut_col = 11
		else:
			cxt.log.write("Peak file {} not in correct format. Should be one of bed, broadpeak, or narrowpeak\n".format(peak_file.fullpath))
			return False
		
		
		bam_files = self.resolve_input('bams', cxt)
		
		frip_output_path = os.path.join(cxt.sample.dest, cxt.sample.name+'.FRiP.txt')
		with open(frip_output_path, 'w') as results_file:
			results_file.write("Sample\tRole\tMapped_Reads\tReads_In_Bed_Regions\tFRiP\n")
			
			for finfo in bam_files:
				cxt.log.write("Computing FRiP for {}\n".format(finfo.fullpath))
				cxt.log.flush()
				
				total_mapped_reads_cmd = 'samtools flagstat "{}" | grep -Po "([0-9]+)(?= \+ [0-9]+ mapped)"'.format(finfo.fullpath)
				cxt.log.write(total_mapped_reads_cmd+'\n')
				cxt.log.flush()
				total_mapped_reads   = float(subprocess.check_output(total_mapped_reads_cmd, shell=True))
				
				reads_in_bed_regions_cmd = 'bedtools coverage -abam "{}" -b "{}" -counts | cut -f {} | paste -sd+ - | bc'.format(finfo.fullpath, peak_file, cut_col)
				cxt.log.write(reads_in_bed_regions_cmd+'\n')
				cxt.log.flush()
				reads_in_bed_regions = float(subprocess.check_output(reads_in_bed_regions_cmd, shell=True))
				frip = reads_in_bed_regions / total_mapped_reads

				results_file.write("%s\t%s\t%d\t%d\t%f\n" % (cxt.sample.name, finfo.cxt.role, total_mapped_reads, reads_in_bed_regions, frip))
		
		return { 'frip': frip_output_path }
	#end run()
#end class FrIPAnalysis