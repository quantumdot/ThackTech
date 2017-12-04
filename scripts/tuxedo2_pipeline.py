#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, AnalysisPipelineCheckpoint, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner
from ThackTech.Pipelines.PipelineModules.HISAT2Align import HISAT2Align




class Tuxedo2PipelineRawSample(PipelineSample):

    def __init__(self, sample, pe_prefix='_R', postfix=""):
        super(Tuxedo2PipelineSample, self).__init__(sample['Basename'], sample['Genome'], sample['Dest'])
        self.set_attribute('PE', ('PE' in sample and sample['PE']))
        self.discover_files(sample['Path'], pe_prefix, postfix)
    #end __init__()
    
    def discover_files(self, path, pe_prefix, postfix):
        files = []
        paths_tried = []
        compressed_extensions = ['.gz', '.bz2', '.zip', '.tar', '.tar.gz']
        if self.get_attribute('PE'):
            base = "{path}/{name}{pre}{read}{post}.{ext}"
            r1 = base.format(path=path, name=self.name, pre=pe_prefix, post=postfix, read=1, ext='fastq')
            r2 = base.format(path=path, name=self.name, pre=pe_prefix, post=postfix, read=2, ext='fastq')
            paths_tried.extend([r1, r2])
            #sys.stderr.write("trying:\n{}\n{}\n".format(r1, r2))
            if os.path.exists(r1) and os.path.exists(r2):
                files.append(FileInfo(r1, FileContext.from_origin('reads'), mate=1))
                files.append(FileInfo(r2, FileContext.from_origin('reads'), mate=2))
            else:
                for ext in compressed_extensions:
                    paths_tried.extend([r1+ext, r2+ext])
                    if os.path.exists(r1+ext) and os.path.exists(r2+ext):
                        files.append(FileInfo(r1+ext, FileContext.from_origin('reads'), mate=1))
                        files.append(FileInfo(r2+ext, FileContext.from_origin('reads'), mate=2))
                        continue
        else:
            base = os.path.join(path, self.name+postfix+'.fastq')
            paths_tried.append(base)
            if os.path.exists(base):
                files.append(FileInfo(base, FileContext.from_origin('reads')))
            else:
                for ext in compressed_extensions:
                    paths_tried.append(base+ext)
                    if os.path.exists(base+ext):
                        files.append(FileInfo(base+ext, FileContext.from_origin('reads')))
                        continue
        if len(files) <= 0:
            raise IOError('Unable to find '+('PE ' if self.get_attribute('PE') else '')+'reads for '+self.name+'\nTried the following paths:\n'+'\n'.join(paths_tried))
        for f in files:
            self.add_file(f)
    #end find_files()
#end Tuxedo2PipelineSample




def main():
	parser = argparse.ArgumentParser()
    parser.add_argument('manifest', help="Manifest file containing sample information in tab separated format. Should contain the following columns (headers must be present): [Path], [Basename], [PE], [Genome], [Dest]")
    #parser.add_argument('--bowtie-version', choices=['1', '2'], default='1', help="Version of bowtie to run")
    
    available_qc_choices = ['ism', 'fqscreen', 'fastqc']
    parser.add_argument('--qc', action='append', default=[], choices=available_qc_choices+['all'], help="Specify which QC pipelines to run after the alignment process completes. SPP is the cross-correlation analysis provided by ccQualityControl/phantompeakqualtools. PBC is the PRC bottleneck/library complexity estimation. ISM computes the distribution of insert size. FPT computes the \"BAM Fingerprint\" using DeepTools bamFingerprint program, and give a good idea of IP strength, especially for TF-like IPs. rpkm will generate a RPKM normalized BigWig from the aligned BAM file. All will run all available QC modules.")
    parser.add_argument('--pe_pre', default='_R', help="Paired-end prefix. String to insert between the file basename and the pair number when searching for read files. If your FASTQ files are names as [reads_R1.fastq, reads_R2.fastq] then use '_R1', or if reads_1.fastq then use '_1'. This option is only used when in paired end mode. default: _R1")
    parser.add_argument('--sample_postfix', default="", help="Postfix to append when looking for read files (ex lane number: '_001')")
    
    parser.add_argument('--unaligned', action='store_true', help='Output reads that fail to align to the reference genome.')
    parser.add_argument('--trim', action='store_true', help="Use trimmomatic to perform adapter clipping.")
    parser.add_argument('--override-dest', action='store', default=None, help="Override the destination read from the sample manifest.")
    
    parser.add_argument('--assembler', action='store', default='stringtie', choices=['stringtie', 'cufflinks'])
    
    performance_group = add_runner_args(parser)
    performance_group.add_argument('--skipalign', action='store_true', help="Skip the alignment process and only run the QC routines. Assumes you have previously aligned files in the proper locations.")

    args, additional_args = parser.parse_known_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices



	#get and run the read alignment pipeline
    pipeline = make_read_alignment_pipeline(args, additional_args)	
	runner = get_configured_runner(args, pipeline)
    runner.run(samples)
	sys.stdout.write("Completed alignment phase for all manifest items!\n")
    sys.stdout.write("=========================================================\n\n")
    sys.stdout.flush()
    
    #get and run the Transcript Merge pipeline
    	#process samples from previous step, and generate new pseudo-sample for transcript mergeing
    
    pipeline = make_transcript_merge_pipeline(args, additional_args)	
	runner = get_configured_runner(args, pipeline)
    runner.run(samples)
	sys.stdout.write("Completed alignment phase for all manifest items!\n")
    sys.stdout.write("=========================================================\n\n")
    sys.stdout.flush()



def make_read_alignment_pipeline(args, additional_args):

    pipeline = AnalysisPipeline('Read Alignment and QC')
    
    if args.trim:
        from ThackTech.Pipelines.PipelineModules import Trimmomatic
        x = Trimmomatic.Trimmomatic(critical=True, processors=args.threads)
        def resolve_trimmomatic_reads(cxt):
            return cxt.sample.find_files(lambda f: f.cxt.role == "reads" )
        x.set_resolver('fastq', resolve_trimmomatic_reads)
        pipeline.append_module(x)
        
        def resolve_bowtie1(cxt):
            if cxt.sample.get_attribute('PE'):
                return cxt.sample.find_files(lambda f: f.cxt.role == "filtered_paired_reads")
            else:
                return cxt.sample.find_files(lambda f: f.cxt.role == "filtered_reads")
    else: 
        def resolve_bowtie1(cxt):
            return cxt.sample.find_files(lambda f: f.cxt.role == "reads" )
    
    
    
    if 'fastqc' in args.qc:
        if args.trim:
            def fastqc_resolver(cxt):
                return [] + resolve_trimmomatic_reads(cxt) + resolve_bowtie1(cxt)
        else:
            def fastqc_resolver(cxt):
                return resolve_bowtie1(cxt)
            
        from ThackTech.Pipelines.PipelineModules import FastQC
        x = FastQC.FastQC(processors=args.threads)
        x.set_resolver('fastqs', fastqc_resolver)
        pipeline.append_module(x)
    
    
    if 'fqscreen' in args.qc:    
        from ThackTech.Pipelines.PipelineModules import FastqScreen
        x = FastqScreen.FastqScreen(processors=args.threads)
        x.set_resolver('fastqs', resolve_bowtie1)
        pipeline.append_module(x)
    
    
    #Align with HISAT2
    #ex: hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
    from ThackTech.Pipelines.PipelineModules import HISAT2Align
    x = HISAT2Align.HISAT2Align(processors=args.threads)
    x.set_parameter('dta', args.assembler)
    x.set_resolver('fastq', resolve_bowtie1)
    pipeline.append_module(x, critical=True)
    
    
    #convert SAM to BAM 
    from ThackTech.Pipelines.PipelineModules import SamToBam
    x = SamToBam.SamToBam(processors=args.threads)
    x.set_resolver('sam', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role == 'sam')[0])
    pipeline.append_module(x, critical=True)
    
    #QC functions:
    if (args.qc is not None) and (len(args.qc) > 0):
    	def qc_bt_bam_resolver(cxt):
            return cxt.sample.find_files(lambda f: f.basename == '{}.bam'.format(cxt.sample.name))[0]
    	
	    if 'ism' in args.qc:
	        from ThackTech.Pipelines.PipelineModules import InsertSizeMetrics
	        x = InsertSizeMetrics.InsertSizeMetrics()
	        x.set_resolver('bam', qc_bt_bam_resolver)
	        pipeline.append_module(x)
	        
    
    #for each sample, assemble transcripts
    #ex: stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf –l ERR188044 ERR188044_chrX.bam
    #dont forget to append the GTF from this step to the mergelist.txt file
    from ThackTech.Pipelines.PipelineModules import StringTie
    x = StringTie.StringTie(processors=args.threads)
    x.set_resolver('alignments', lambda cxt: return None)

	return pipeline
#end make_read_alignment_pipeline()
    
   
    
    
def make_transcript_merge_pipeline(args):
	pipeline = AnalysisPipeline('Transcript Merge') 
    #########################
    #
    # Start New Pipeline
    #
    #########################
     
    # Merge Transcript assemblies for all samples
    #probably requires only one node
    
    #optional: Examine how the transcripts compare with the reference annotation
    #ex: gffcompare –r chrX_data/genes/chrX.gtf –G –o merged stringtie_merged.gtf
    
#end make_transcript_merge_pipeline()

   
    
def make_transcript_quant_pipeline(args):   
    #########################
    #
    # Start New Pipeline
    #
    #########################
    
    #Estimate transcript abundances and create table counts for Ballgown:
    #stringtie –e –B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
    
    
    
#end make_transcript_quant_pipeline()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

