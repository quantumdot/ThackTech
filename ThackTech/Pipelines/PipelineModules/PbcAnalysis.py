import os
import subprocess
import shlex
from ThackTech import filetools, Common
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class PbcAnalysis(PipelineModule):
	
	def __init__(self, **kwargs):
		super(PbcAnalysis, self).__init__('PBC', 'Cross-correlation analysis using SPP', **kwargs)
		self._name_resolver('bam')
	#end __init__()
	
	def run(self, cxt):
		cxt.log.write('\t-> Running PBC QC...\n')
		cxt.log.flush()
		#compute PBC
		bam = self.resolve_input('bam', cxt.sample)[0].fullpath
		results = self.run_filter_qc(bam, cxt.sample.get_attribute('PE'), "-q 30")
		pbc_dir = os.path.join(cxt.sample.dest, 'pbc')
		filetools.ensure_dir(pbc_dir)
		self._run_subprocess('mv -f '+cxt.sample.name+'.filt.* '+cxt.sample.name+'.dup.qc '+pbc_dir, shell=True)
		self._run_subprocess(['rm', '-f', os.path.join(cxt.sample.dest), ('tmp.%s.filt.srt.bam' % (cxt.sample.name,))])
		cxt.log.write('\t-> Completed PBC QC...\n')
		cxt.log.flush()
		
		output_files = []
		output_files.append(FileInfo(os.path.join(pbc_dir, '{}.filt.srt.bam'.format(cxt.sample.name)),
							FileContext.from_module_context(cxt, 'filtered_bam')))
		
		output_files.append(FileInfo(os.path.join(pbc_dir, '{}.filt.srt.nodup.bam'.format(cxt.sample.name)),
							FileContext.from_module_context(cxt, 'filtered_deduplicated_bam')))
		
		output_files[1].companions.append(FileInfo(os.path.join(pbc_dir, '{}.filt.srt.nodup.bai'.format(cxt.sample.name)),
										  FileContext.from_module_context(cxt, 'filtered_deduplicated_bam_idx')))
		
		output_files.append(FileInfo(os.path.join(pbc_dir, '{}.filt.srt.nodup.flagstat.qc'.format(cxt.sample.name)),
							FileContext.from_module_context(cxt, 'filtered_deduplicated_flagstat')))
		
		output_files.append(FileInfo(os.path.join(pbc_dir, '{}.filt.srt.nodup.pbc.qc'.format(cxt.sample.name)),
							FileContext.from_module_context(cxt, 'duplicate_qc')))
		
		output_files.append(FileInfo(os.path.join(pbc_dir, '{}.dup.qc'.format(cxt.sample.name)),
							FileContext.from_module_context(cxt, 'pbc_qc')))

		return output_files
	#end run()
	
	def run_filter_qc(self, input_bam=None, paired_end=None, samtools_params=None):
		if input_bam is None:
			raise Exception('input_bam is required!')
	
		if paired_end is None:
			raise Exception('paired_end is required!')
	
		
	
		# The following line(s) download your file inputs to the local file system
		# using variable names for the filenames.
		input_bam = os.path.abspath(input_bam)
		raw_bam_dir = os.path.dirname(input_bam)
		raw_bam_filename = os.path.basename(input_bam)
		raw_bam_basename = os.path.splitext(raw_bam_filename)[0]
		origional_dir = os.getcwd()
		os.chdir(raw_bam_dir)
		
		filt_bam_prefix = raw_bam_basename + ".filt.srt" 
		filt_bam_filename = filt_bam_prefix + ".bam"
		if paired_end:
			# =============================
			# Remove  unmapped, mate unmapped
			# not primary alignment, reads failing platform
			# Remove low MAPQ reads
			# Only keep properly paired reads
			# Obtain name sorted BAM file
			# ==================
			tmp_filt_bam_prefix = "tmp.%s" %(filt_bam_prefix) #was tmp.prefix.nmsrt
			tmp_filt_bam_filename = tmp_filt_bam_prefix + ".bam"
			out, err = Common.run_pipe([
				#filter:  -F 1804 FlAG bits to exclude; -f 2 FLAG bits to reqire; -q 30 exclude MAPQ < 30; -u uncompressed output
				#exclude FLAG 1804: unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
				#require FLAG 2: properly aligned
				"samtools view -F 1804 -f 2 %s -u %s" % (samtools_params, raw_bam_filename),
				#sort:  -n sort by name; - take input from stdin; out to specified filename
				"samtools sort -n - %s" % (tmp_filt_bam_prefix)])  # Will produce name sorted BAM
			if err:
				raise RuntimeError("samtools error: %s" %(err))
			
			
			# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
			# Obtain position sorted BAM
			#print subprocess.check_output('ls -l', shell=True)
			out,err = Common.run_pipe([
				#fill in mate coordinates, ISIZE and mate-related flags
				#fixmate requires name-sorted alignment; -r removes secondary and unmapped (redundant here because already done above?)
				#- send output to stdout
				"samtools fixmate -r %s -" %(tmp_filt_bam_filename),
				#repeat filtering after mate repair
				"samtools view -F 1804 -f 2 -u -",
				#produce the coordinate-sorted BAM
				"samtools sort - %s" %(filt_bam_prefix)])
			#print subprocess.check_output('ls -l', shell=True)
		else: #single-end data
			# =============================
			# Remove unmapped, mate unmapped
			# not primary alignment, reads failing platform
			# Remove low MAPQ reads
			# Obtain name sorted BAM file
			# ==================  
			with open(filt_bam_filename, 'w') as fh:
				subprocess.check_call(['samtools', 'view', '-F', '1804', samtools_params, '-b', raw_bam_filename], stdout=fh)
		
		# ========================
		# Mark duplicates
		# ======================
		tmp_filt_bam_filename = raw_bam_basename + ".dupmark.bam"
		dup_file_qc_filename = raw_bam_basename + ".dup.qc"
		subprocess.check_call(shlex.split(
		 "picard-tools MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s \
		  VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
		  %(filt_bam_filename, tmp_filt_bam_filename, dup_file_qc_filename)))
		os.rename(tmp_filt_bam_filename,filt_bam_filename)
		
		if paired_end:
			final_bam_prefix = raw_bam_basename + ".filt.srt.nodup"
		else:
			final_bam_prefix = raw_bam_basename + ".filt.nodup.srt"
		final_bam_filename = final_bam_prefix + ".bam" # To be stored
		final_bam_index_filename = final_bam_prefix + ".bai" # To be stored
		final_bam_file_mapstats_filename = final_bam_prefix + ".flagstat.qc" # QC file
		
		if paired_end:
			# ============================
			# Remove duplicates
			# Index final position sorted BAM
			# Create final name sorted BAM
			# ============================
			with open(final_bam_filename, 'w') as fh:
				subprocess.check_call(['samtools', 'view', '-F', '1804', '-f2', '-b', filt_bam_filename], stdout=fh)
			#namesorting is needed for bam->bedPE, so moved to xcor
			#final_nmsrt_bam_prefix = raw_bam_basename + ".filt.nmsrt.nodup"
			#final_nmsrt_bam_filename = final_nmsrt_bam_prefix + ".bam"
			#subprocess.check_call(shlex.split("samtools sort -n %s %s" %(final_bam_filename, final_nmsrt_bam_prefix)))
		else:
			# ============================
			# Remove duplicates
			# Index final position sorted BAM
			# ============================
			with open(final_bam_filename, 'w') as fh:
				subprocess.check_call(['samtools', 'view', '-F', '1804', '-b', filt_bam_filename], stdout=fh)
		# Index final bam file
		subprocess.check_call(['samtools', 'index', final_bam_filename, final_bam_index_filename])
		# Generate mapping statistics
		with open(final_bam_file_mapstats_filename, 'w') as fh:
			subprocess.check_call(['samtools', 'flagstat', final_bam_filename], stdout=fh)
		
		# =============================
		# Compute library complexity
		# =============================
		# Sort by name
		# convert to bedPE and obtain fragment coordinates
		# sort by position and strand
		# Obtain unique count statistics
		pbc_file_qc_filename = final_bam_prefix + ".pbc.qc"
		# PBC File output
		# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
		if paired_end:
			steps = [
				"samtools sort -no %s -" %(filt_bam_filename),
				"bamToBed -bedpe -i stdin",
				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}'"""
			]
		else:
			steps = [
				"bamToBed -i %s" %(filt_bam_filename), #for some reason 'bedtools bamtobed' does not work but bamToBed does
				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}'"""
			]
		# these st
		steps.extend([
			"grep -v 'chrM'", #TODO this should be implemented as an explicit list of allowable names, so that mapping can be done to a complete reference
			"sort",
			"uniq -c",
			r"""awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'"""
		 ])
		out,err = Common.run_pipe(steps,pbc_file_qc_filename)
		if err:
			raise Exception("PBC file error: %s" %(err))
		
		headers = 'TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair'
		p = subprocess.Popen('echo "'+headers+'" | cat - "'+pbc_file_qc_filename+'" > /tmp/out && mv /tmp/out "'+pbc_file_qc_filename+'"', shell=True)
		p.communicate()
		
		# Return links to the output files
		output = {
			"filtered_bam": final_bam_filename,
			"filtered_bam_index": final_bam_index_filename,
			"filtered_mapstats": final_bam_file_mapstats_filename,
			"dup_file_qc": dup_file_qc_filename,
			"pbc_file_qc": pbc_file_qc_filename,
			"paired_end": paired_end
		}
		#output.update({'output_JSON': output.copy()})
		os.chdir(origional_dir)
		#print "Exiting with output: %s" %(output)
		return output
	#end run_filter_qc()
		
		
		
#end PbcAnalysis




	