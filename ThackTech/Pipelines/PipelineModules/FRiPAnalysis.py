import os
#import subprocess
#from ThackTech import filetools
from ThackTech.Pipelines import PipelineModule
from deeptools import countReadsPerBin
import pysam

class FRiPAnalysis(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='FRiP', short_description='Fraction of Reads in Peaks')
		super_args.update(**kwargs)
		super(FRiPAnalysis, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bed')
		self._name_resolver('bams')
	#end __declare_resolvers()
	
	
	def run(self, cxt):
		peak_file = self.resolve_input('bed', cxt)
		bam_files = self.resolve_input('bams', cxt)
		
		with open(peak_file.fullpath, 'r') as bed_file:
			cr = countReadsPerBin.CountReadsPerBin([f.fullpath for f in bam_files],
			                                        bedFile=bed_file,
			                                        numberOfProcessors=self.processors)
			reads_at_peaks = cr.run()
			sum_reads_in_peaks = reads_at_peaks.sum(axis=0)
			
		frip_output_path = os.path.join(cxt.sample.dest, cxt.sample.name+'.FRiP.txt')
		with open(frip_output_path, 'w') as results_file:	
			results_file.write("Sample\tRole\tMapped_Reads\tReads_In_Bed_Regions\tFRiP\n")
			
			for i in range(len(bam_files)):
				finfo = bam_files[i]
				b = pysam.AlignmentFile(f.fullpath)
				results_file.write("{name}\t{role}\t{mapped}\t{overlap}\t{frip}\n".format(name=cxt.sample.name, 
																						  role=finfo.cxt.role,
																						  mapped=b.mapped, 
																						  overlap=sum_reads_in_peaks[i], 
																						  frip=(float(sum_reads_in_peaks[i]) / b.mapped)))
		return { 'frip': frip_output_path }
	#end run2()
	
# 	def run(self, cxt):
# 		return self.run2(cxt)
# 		 
# 		
# 		peak_file = self.resolve_input('bed', cxt)
# 		if not peak_file.isfile:
# 			raise IOError('Unable to locate peak file {}!'.format(peak_file.fullpath))
# 		peak_type = peak_file.ext.lower()
# 		if peak_type == '.bed':
# 			cut_col = 6
# 		elif peak_type == '.broadpeak':
# 			cut_col = 10
# 		elif peak_type == '.narrowpeak':
# 			cut_col = 11
# 		else:
# 			cxt.log.write("Peak file {} not in correct format. Should be one of bed, broadpeak, or narrowpeak\n".format(peak_file.fullpath))
# 			return False
# 		
# 		
# 		bam_files = self.resolve_input('bams', cxt)
# 		
# 		frip_output_path = os.path.join(cxt.sample.dest, cxt.sample.name+'.FRiP.txt')
# 		with open(frip_output_path, 'w') as results_file:
# 			results_file.write("Sample\tRole\tMapped_Reads\tReads_In_Bed_Regions\tFRiP\n")
# 			
# 			for finfo in bam_files:
# 				cxt.log.write("Computing FRiP for {}\n".format(finfo.fullpath))
# 				cxt.log.flush()
# 				
# 				total_mapped_reads_cmd = 'samtools flagstat "{}" | grep -Po "([0-9]+)(?= \+ [0-9]+ mapped)"'.format(finfo.fullpath)
# 				cxt.log.write(total_mapped_reads_cmd+'\n')
# 				cxt.log.flush()
# 				total_mapped_reads   = float(subprocess.check_output(total_mapped_reads_cmd, shell=True))
# 				
# 				reads_in_bed_regions_cmd = 'bedtools coverage -a "{}" -b "{}" -counts | cut -f {} | paste -sd+ - | bc'.format(peak_file.fullpath, finfo.fullpath, cut_col)
# 				cxt.log.write(reads_in_bed_regions_cmd+'\n')
# 				cxt.log.flush()
# 				reads_in_bed_regions = float(subprocess.check_output(reads_in_bed_regions_cmd, shell=True))
# 				frip = reads_in_bed_regions / total_mapped_reads
# 
# 				results_file.write("%s\t%s\t%d\t%d\t%f\n" % (cxt.sample.name, finfo.cxt.role, total_mapped_reads, reads_in_bed_regions, frip))
# 		
# 		return { 'frip': frip_output_path }
# 	#end run()
#end class FrIPAnalysis