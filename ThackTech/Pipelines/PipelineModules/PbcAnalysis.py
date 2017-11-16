import os
import re
import shlex
import subprocess
from pprint import pprint, pformat
from ThackTech import filetools, Common
from ThackTech.Pipelines import PipelineModule
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext


class PbcAnalysis(PipelineModule):
	
	def __init__(self, **kwargs):
		super_args = dict(name='PBC', short_description='PCR Bottleneck Analysis')
		super_args.update(**kwargs)
		super(PbcAnalysis, self).__init__(**super_args)
	#end __init__()
	
	def _declare_parameters(self):
		pass
	#end __declare_parameters()
	
	def _declare_resolvers(self):
		self._name_resolver('bam')
	#end __declare_resolvers()
	
	def run(self, cxt):
		cxt.log.write('\t-> Running PBC QC...\n')
		cxt.log.flush()
		#compute PBC
		bam = self.resolve_input('bam', cxt)
		
		#set up an isolated directory to work in
		origional_dir = os.getcwd()
		tmp_pbc_dir = os.path.join(cxt.sample.dest, 'pbc_'+cxt.sample.name)
		filetools.ensure_dir(tmp_pbc_dir)
		os.chdir(tmp_pbc_dir)
		
		#Run the PBC analysis
		results = self.run_filter_qc2(cxt, bam, "-q 30")
		
		#restore the origional working directory
		os.chdir(origional_dir)
		
		final_pbc_dir = os.path.join(cxt.sample.dest, 'pbc')
		filetools.ensure_dir(final_pbc_dir)
		for rf in results:
			rf.move(final_pbc_dir)
		
		self._run_subprocess(['rm', '-rf', tmp_pbc_dir])
		cxt.log.write('\t-> Completed PBC QC...\n')
		cxt.log.flush()

		return results
	#end run()
	

	
	
	
	def run_filter_qc2(self, cxt, input_bam, samtools_params=None):
		''' Run the PBC analysis
		
			Parameters:
				cxt: ModuleRunContext
				input_bam: FileInfo representing the bam to operate on
				
		'''
		
		#Set this to indicate PE or not
		paired_end = cxt.sample.get_attribute('PE')

		# Generate initial mapping statistics
		cxt.log.write("Generating initial mapping statistics\n")
		raw_bam_file_mapstats_filename = input_bam.basename_with_ext('flagstat.qc')
		with open(raw_bam_file_mapstats_filename, 'w') as fh:
			flagstat_command = "samtools flagstat {}".format(input_bam.fullpath)
			cxt.log.write(flagstat_command+"\n")
			subprocess.check_call(shlex.split(flagstat_command), stdout=fh)

		filt_bam_prefix = input_bam.basename_with_ext("filt.srt")
		filt_bam_filename = filt_bam_prefix + ".bam"
		if paired_end:
			# =============================
			# Remove  unmapped, mate unmapped
			# not primary alignment, reads failing platform
			# Remove low MAPQ reads
			# Only keep properly paired reads
			# Obtain name sorted BAM file
			# =============================
			tmp_filt_bam_prefix = "tmp.{}".format(filt_bam_prefix)  # was tmp.prefix.nmsrt
			tmp_filt_bam_filename = tmp_filt_bam_prefix + ".bam"
			out, err = Common.run_pipe([
				# filter: -F 1804 FlAG bits to exclude; -f 2 FLAG bits to require;
				# -q 30 exclude MAPQ < 30; -u uncompressed output
				# exclude FLAG 1804: unmapped, next segment unmapped, secondary
				# alignments, not passing platform q, PCR or optical duplicates
				# require FLAG 2: properly aligned
				"samtools view -F 1804 -f 2 {stparams} -u {bam}".format(stparams=samtools_params, bam=input_bam.fullpath),
				
				# sort:  -n sort by name; - take input from stdin;
				# out to specified filename
				# Will produce name sorted BAM
				"samtools sort -n - {}".format(tmp_filt_bam_prefix)])
			if err:
				cxt.log.write("samtools error: {}\n".format(err))
				
				
			# Remove orphan reads (pair was removed)and read pairs mapping to different chromosomes
			# Obtain position sorted BAM
			out, err = Common.run_pipe([
				# fill in mate coordinates, ISIZE and mate-related flags
				# fixmate requires name-sorted alignment;
				# - send output to stdout
				"samtools fixmate {} -".format(tmp_filt_bam_filename),
				
				# repeat filtering after mate repair
				"samtools view -F 1804 -f 2 -u -",
				
				# produce the coordinate-sorted BAM
				"samtools sort - {}".format(filt_bam_prefix)])
			
		else:  # single-end data
			# =============================
			# Remove unmapped, mate unmapped
			# not primary alignment, reads failing platform
			# Remove low MAPQ reads
			# Obtain name sorted BAM file
			# ==================
			with open(filt_bam_filename, 'w') as fh:
				samtools_filter_command = "samtools view -F 1804 {stparams} -b {bam}".format(stparams=samtools_params, bam=input_bam.fullpath)
				cxt.log.write(samtools_filter_command+"\n")
				subprocess.check_call(shlex.split(samtools_filter_command), stdout=fh)

		# ========================
		# Mark duplicates
		# ======================
		tmp_filt_bam_filename = input_bam.basename_with_ext("dupmark.bam")
		dup_file_qc_filename = input_bam.basename_with_ext("dup.qc")
		picard_string = ' '.join([
			"picard-tools MarkDuplicates",
			"INPUT=%s" % (filt_bam_filename),
			"OUTPUT=%s" % (tmp_filt_bam_filename),
			"METRICS_FILE=%s" % (dup_file_qc_filename),
			"VALIDATION_STRINGENCY=LENIENT",
			"ASSUME_SORTED=true",
			"REMOVE_DUPLICATES=false"
			])
		cxt.log.write(picard_string+"\n")
		subprocess.check_output(shlex.split(picard_string))
		os.rename(tmp_filt_bam_filename, filt_bam_filename)

		if paired_end:
			final_bam_prefix = input_bam.basename_with_ext("filt.srt.nodup")
		else:
			final_bam_prefix = input_bam.basename_with_ext(".filt.nodup.srt")
			
		final_bam_filename = final_bam_prefix + ".bam"  # To be stored
		final_bam_index_filename = final_bam_filename + ".bai"  # To be stored
		# QC file
		final_bam_file_mapstats_filename = final_bam_prefix + ".flagstat.qc"

		if paired_end:
			samtools_dedupe_command = "samtools view -F 1804 -f2 -b %s" % (filt_bam_filename)
		else:
			samtools_dedupe_command = "samtools view -F 1804 -b %s" % (filt_bam_filename)

		# ============================
		# Remove duplicates
		# Index final position sorted BAM
		# ============================
		with open(final_bam_filename, 'w') as fh:
			cxt.log.write(samtools_dedupe_command+"\n")
			subprocess.check_call(shlex.split(samtools_dedupe_command), stdout=fh)
			
		# Index final bam file
		samtools_index_command = "samtools index {} {}".format(final_bam_filename, final_bam_index_filename)
		cxt.log.write(samtools_index_command+"\n")
		subprocess.check_output(shlex.split(samtools_index_command))

		# Generate mapping statistics
		with open(final_bam_file_mapstats_filename, 'w') as fh:
			flagstat_command = "samtools flagstat {}".format(final_bam_filename)
			cxt.log.write(flagstat_command+"\n")
			subprocess.check_call(shlex.split(flagstat_command), stdout=fh)

		# =============================
		# Compute library complexity
		# =============================
		# Sort by name
		# convert to bedPE and obtain fragment coordinates
		# sort by position and strand
		# Obtain unique count statistics
		pbc_file_qc_filename = final_bam_prefix + ".pbc.qc"
		# PBC File output
		# TotalReadPairs [tab]
		# DistinctReadPairs [tab]
		# OneReadPair [tab]
		# TwoReadPairs [tab]
		# NRF=Distinct/Total [tab]
		# PBC1=OnePair/Distinct [tab]
		# PBC2=OnePair/TwoPair
		if paired_end:
			steps = [
				"samtools sort -no {} -".format(filt_bam_filename),
				"bamToBed -bedpe -i stdin",
				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}'"""]
		else:
			steps = [
				"bamToBed -i {}".format(filt_bam_filename),
				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}'"""]
		
		steps.extend([
			"grep -v 'chrM'",
			"sort",
			"uniq -c",
			r"""awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'"""
			])
		out, err = Common.run_pipe(steps, pbc_file_qc_filename)
		if err:
			cxt.log.write("PBC file error: {}\n".format(err))

		
		cxt.log.write("Calculating QC metrics\n")
		dup_qc = dup_parse(dup_file_qc_filename)
		pbc_qc = pbc_parse(pbc_file_qc_filename)
		initial_mapstats_qc = flagstat_parse(raw_bam_file_mapstats_filename)
		final_mapstats_qc = flagstat_parse(final_bam_file_mapstats_filename)
		
		if paired_end:
			useable_fragments = final_mapstats_qc.get('in_total')[0]/2
		else:
			useable_fragments = final_mapstats_qc.get('in_total')[0]
			
		cxt.log.write("initial_mapstats_qc: {}\n".format(initial_mapstats_qc)),
		cxt.log.write("final_mapstats_qc: {}\n".format(final_mapstats_qc)),
		cxt.log.write("dup_qc: {}\n".format(dup_qc))
		cxt.log.write("pbc_qc: {}\n".format(pbc_qc))

		# Return links to the output files and values.
		out_files = []
		out_files.append(FileInfo(final_bam_filename, FileContext.from_module_context(cxt, "filtered_bam")))
		out_files[0].companions.append(FileInfo(final_bam_index_filename, FileContext.from_module_context(cxt, "filtered_bam_index")))
		out_files.append(FileInfo(final_bam_file_mapstats_filename, FileContext.from_module_context(cxt, "filtered_mapstats")))
		out_files.append(FileInfo(dup_file_qc_filename, FileContext.from_module_context(cxt, "dup_file_qc")))
		out_files.append(FileInfo(pbc_file_qc_filename, FileContext.from_module_context(cxt, "pbc_file_qc")))
		
		
		data = {
			"paired_end": paired_end,
			"n_reads_input": str(initial_mapstats_qc.get('in_total')[0]),
			"picard_read_pairs_examined": str(dup_qc.get('read_pairs_examined')),
			"picard_unpaired_reads_examined": str(dup_qc.get('unpaired_reads_examined')),
			"picard_read_pair_duplicates": str(dup_qc.get('read_pair_duplicates')),
			"picard_unpaired_read_duplicates": str(dup_qc.get('unpaired_read_duplicates')),
			"useable_fragments": str(useable_fragments),
			"NRF": str(pbc_qc.get('NRF')),
			"PBC1": str(pbc_qc.get('PBC1')),
			"PBC2": str(pbc_qc.get('PBC2')),
			"duplicate_fraction": str(dup_qc.get('percent_duplication'))
		}
		cxt.log.write("Exiting with output:\n{}\n\n".format(pformat(data)))
		
		return out_files
	#end run_filter_qc2()
	
	
	
	
	
	
	
# 	def run_filter_qc(self, cxt, input_bam=None, paired_end=None, samtools_params=None):
# 		if input_bam is None:
# 			raise Exception('input_bam is required!')
# 	
# 		if paired_end is None:
# 			raise Exception('paired_end is required!')
# 	
# 		#filter_exclude = 1548 #read unmapped, mate unmapped, read fails platform/vendor quality checksm, read is PCR or optical duplicate
# 		filter_exclude = 1804 #read unmapped, mate unmapped, read fails platform/vendor quality checksm, read is PCR or optical duplicate, not primary alignment
# 
# 	
# 		# The following line(s) download your file inputs to the local file system
# 		# using variable names for the filenames.
# 		input_bam = os.path.abspath(input_bam)
# 		raw_bam_dir = os.path.dirname(input_bam)
# 		raw_bam_filename = os.path.basename(input_bam)
# 		raw_bam_basename = os.path.splitext(raw_bam_filename)[0]
# 		origional_dir = os.getcwd()
# 		os.chdir(raw_bam_dir)
# 		
# 		filt_bam_prefix = raw_bam_basename + ".filt.srt" 
# 		filt_bam_filename = filt_bam_prefix + ".bam"
# 		if paired_end:
# 			# =============================
# 			# Remove  unmapped, mate unmapped
# 			# not primary alignment, reads failing platform
# 			# Remove low MAPQ reads
# 			# Only keep properly paired reads
# 			# Obtain name sorted BAM file
# 			# ==================
# 			tmp_filt_bam_prefix = "tmp.%s" %(filt_bam_prefix) #was tmp.prefix.nmsrt
# 			tmp_filt_bam_filename = tmp_filt_bam_prefix + ".bam"
# 			out, err = Common.run_pipe([
# 				#filter:  -F 1804 FlAG bits to exclude; -f 2 FLAG bits to reqire; -q 30 exclude MAPQ < 30; -u uncompressed output
# 				#exclude FLAG 1804: unmapped, next segment unmapped, secondary alignments, not passing platform q, PCR or optical duplicates
# 				#require FLAG 2: properly aligned
# 				"samtools view -F %d -f 2 %s -u %s" % (filter_exclude, samtools_params, raw_bam_filename),
# 				#sort:  -n sort by name; - take input from stdin; out to specified filename
# 				"samtools sort -n - %s" % (tmp_filt_bam_prefix)], stderr=cxt.log)  # Will produce name sorted BAM
# 			if err:
# 				raise RuntimeError("samtools error: %s" %(err))
# 			
# 			
# 			# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
# 			# Obtain position sorted BAM
# 			#print subprocess.check_output('ls -l', shell=True)
# 			out,err = Common.run_pipe([
# 				#fill in mate coordinates, ISIZE and mate-related flags
# 				#fixmate requires name-sorted alignment; -r removes secondary and unmapped (redundant here because already done above?)
# 				#- send output to stdout
# 				"samtools fixmate %s -" %(tmp_filt_bam_filename),
# 				#repeat filtering after mate repair
# 				"samtools view -F %d -f 2 -u -" % (filter_exclude,),
# 				#produce the coordinate-sorted BAM
# 				"samtools sort - %s" %(filt_bam_prefix)], stderr=cxt.log)
# 			#print subprocess.check_output('ls -l', shell=True)
# 		else: #single-end data
# 			# =============================
# 			# Remove unmapped, mate unmapped
# 			# not primary alignment, reads failing platform
# 			# Remove low MAPQ reads
# 			# Obtain name sorted BAM file
# 			# ==================  
# 			with open(filt_bam_filename, 'w') as fh:
# 				subprocess.check_call(['samtools', 'view', '-F', str(filter_exclude), samtools_params, '-b', raw_bam_filename], stdout=fh, stderr=cxt.log)
# 				
# 		subprocess.check_call(['samtools', 'index', filt_bam_filename, filt_bam_prefix+'.bai'], stderr=cxt.log)
# 		
# 		# ========================
# 		# Mark duplicates
# 		# ======================
# 		tmp_filt_bam_filename = raw_bam_basename + ".dupmark.bam"
# 		dup_file_qc_filename = raw_bam_basename + ".dup.qc"
# 		subprocess.check_call(shlex.split(
# 		 "picard-tools MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s \
# 		  VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
# 		  %(filt_bam_filename, tmp_filt_bam_filename, dup_file_qc_filename)), stderr=cxt.log)
# 		os.rename(tmp_filt_bam_filename,filt_bam_filename)
# 		
# 		final_bam_prefix = raw_bam_basename + ".filt.srt.nodup"
# 		final_bam_filename = final_bam_prefix + ".bam" # To be stored
# 		final_bam_index_filename = final_bam_prefix + ".bai" # To be stored
# 		final_bam_file_mapstats_filename = final_bam_prefix + ".flagstat.qc" # QC file
# 		
# 		if paired_end:
# 			# ============================
# 			# Remove duplicates
# 			# Index final position sorted BAM
# 			# Create final name sorted BAM
# 			# ============================
# 			with open(final_bam_filename, 'w') as fh:
# 				subprocess.check_call(['samtools', 'view', '-F', str(filter_exclude), '-f2', '-b', filt_bam_filename], stdout=fh, stderr=cxt.log)
# 			#namesorting is needed for bam->bedPE, so moved to xcor
# 			#final_nmsrt_bam_prefix = raw_bam_basename + ".filt.nmsrt.nodup"
# 			#final_nmsrt_bam_filename = final_nmsrt_bam_prefix + ".bam"
# 			#subprocess.check_call(shlex.split("samtools sort -n %s %s" %(final_bam_filename, final_nmsrt_bam_prefix)))
# 		else:
# 			# ============================
# 			# Remove duplicates
# 			# Index final position sorted BAM
# 			# ============================
# 			with open(final_bam_filename, 'w') as fh:
# 				subprocess.check_call(['samtools', 'view', '-F', str(filter_exclude), '-b', filt_bam_filename], stdout=fh, stderr=cxt.log)
# 		# Index final bam file
# 		subprocess.check_call(['samtools', 'index', final_bam_filename, final_bam_index_filename], stderr=cxt.log)
# 		# Generate mapping statistics
# 		with open(final_bam_file_mapstats_filename, 'w') as fh:
# 			subprocess.check_call(['samtools', 'flagstat', final_bam_filename], stdout=fh, stderr=cxt.log)
# 		
# 		# =============================
# 		# Compute library complexity
# 		# =============================
# 		# Sort by name
# 		# convert to bedPE and obtain fragment coordinates
# 		# sort by position and strand
# 		# Obtain unique count statistics
# 		pbc_file_qc_filename = final_bam_prefix + ".pbc.qc"
# 		# PBC File output
# 		# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
# 		if paired_end:
# 			steps = [
# 				"samtools sort -no %s %s" % (filt_bam_filename, 'tmp_name_srt.'+filt_bam_filename),
# 				"bamToBed -bedpe -i stdin",
# 				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}'"""
# 			]
# 		else:
# 			steps = [
# 				"bamToBed -i %s" % (filt_bam_filename), #for some reason 'bedtools bamtobed' does not work but bamToBed does
# 				r"""awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}'"""
# 			]
# 		# these st
# 		steps.extend([
# 			"grep -v 'chrM'", #TODO this should be implemented as an explicit list of allowable names, so that mapping can be done to a complete reference
# 			"sort",
# 			"uniq -c",
# 			r"""awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'"""
# 		 ])
# 		out,err = Common.run_pipe(steps, outfile=pbc_file_qc_filename, stderr=cxt.log)
# 		if err:
# 			raise Exception("PBC file error: %s" %(err))
# 		
# 		headers = 'TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair\n'
# 		filetools.prepend_file(pbc_file_qc_filename, headers)
# 		
# 		# Return links to the output files
# 		output = {
# 			"filtered_bam": final_bam_filename,
# 			"filtered_bam_index": final_bam_index_filename,
# 			"filtered_mapstats": final_bam_file_mapstats_filename,
# 			"dup_file_qc": dup_file_qc_filename,
# 			"pbc_file_qc": pbc_file_qc_filename,
# 			"paired_end": paired_end
# 		}
# 		#output.update({'output_JSON': output.copy()})
# 		os.chdir(origional_dir)
# 		#print "Exiting with output: %s" %(output)
# 		return output
# 	#end run_filter_qc()
		
#end PbcAnalysis



def dup_parse(fname):
	with open(fname, 'r') as dup_file:
		if not dup_file:
			return None
		
		lines = iter(dup_file.read().splitlines())
		
		for line in lines:
			if line.startswith('## METRICS CLASS'):
				headers = lines.next().rstrip('\n').lower()
				metrics = lines.next().rstrip('\n')
				break
		
		headers = headers.split('\t')
		metrics = metrics.split('\t')
		headers.pop(0)
		metrics.pop(0)
		
		dup_qc = dict(zip(headers, metrics))
	return dup_qc
#end dup_parse()

def pbc_parse(fname):
	with open(fname, 'r') as pbc_file:
		if not pbc_file:
			return None

		lines = pbc_file.read().splitlines()
		line = lines[0].rstrip('\n')
		# PBC File output:
		#   TotalReadPairs <tab>
		#   DistinctReadPairs <tab>
		#   OneReadPair <tab>
		#   TwoReadPairs <tab>
		#   NRF=Distinct/Total <tab>
		#   PBC1=OnePair/Distinct <tab>
		#   PBC2=OnePair/TwoPair

		headers = ['TotalReadPairs',	
                   'DistinctReadPairs',
                   'OneReadPair',
                   'TwoReadPairs',
                   'NRF',
                   'PBC1',
                   'PBC2']
		metrics = line.split('\t')

		pbc_qc = dict(zip(headers, metrics))
	return pbc_qc
#end pbc_parse()

def flagstat_parse(fname):
	with open(fname, 'r') as flagstat_file:
		if not flagstat_file:
			return None
		flagstat_lines = flagstat_file.read().splitlines()

	qc_dict = {
        # values are regular expressions,
        # will be replaced with scores [hiq, lowq]
        'in_total': 'in total',
        'duplicates': 'duplicates',
        'mapped': 'mapped',
        'paired_in_sequencing': 'paired in sequencing',
        'read1': 'read1',
        'read2': 'read2',
        'properly_paired': 'properly paired',
        'with_self_mate_mapped': 'with itself and mate mapped',
        'singletons': 'singletons',
        # i.e. at the end of the line
        'mate_mapped_different_chr': 'with mate mapped to a different chr$',
        # RE so must escape
        'mate_mapped_different_chr_hiQ': 'with mate mapped to a different chr \(mapQ>=5\)'
    }

	for (qc_key, qc_pattern) in qc_dict.items():
		qc_metrics = next(re.split(qc_pattern, line)
                          for line in flagstat_lines
                          if re.search(qc_pattern, line))
		(hiq, lowq) = qc_metrics[0].split(' + ')
		qc_dict[qc_key] = [int(hiq.rstrip()), int(lowq.rstrip())]

	return qc_dict
#end flagstat_parse()



	