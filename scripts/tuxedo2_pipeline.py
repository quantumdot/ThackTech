#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, AnalysisPipelineCheckpoint, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner

#sample manifest should be in the following TAB separated format (with headers):
#Path    Basename    PE
#/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true

#BEFORE RUNNING SCRIPT WHEN USING LMOD
#module load java/1.8.0_121 samtools/0.1.19 intel/17.0.2 python/2.7.12 bedtools2/2.25.0 R-Project/3.3.3 HISAT2/2.1.0

class Tuxedo2PipelineRawSample(PipelineSample):

    def __init__(self, sample, genome, dest, pe_prefix='_R', postfix=""):
        super(Tuxedo2PipelineRawSample, self).__init__(sample['Basename'], genome, dest)
        self.set_attribute('PE', ('PE' in sample and sample['PE']))
        
        try:
            self.read_file_manifest(self.default_file_manifest_location)
        except:
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

class Tuxedo2PipelineMergeSample(PipelineSample):
    def __init__(self, sample):
        super(Tuxedo2PipelineMergeSample, self).__init__(sample['Basename'], sample['Genome'], sample['Dest'])
    #end __init__()
#end Tuxedo2PipelineMergeSample


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', help="Manifest file containing sample information in tab separated format. Should contain the following columns (headers must be present): [Path], [Basename], [PE], [Genome], [Dest]")
    parser.add_argument('--dest', action='store', required=True, help="Path to destination for results.")
    parser.add_argument('--genome', action='store', required=True, help="Reference genome to use.")
    
    available_qc_choices = ['ism', 'fqscreen', 'fastqc']
    parser.add_argument('--qc', action='append', default=[], choices=available_qc_choices+['all'], help="Specify which QC pipelines to run after the alignment process completes.")
    parser.add_argument('--pe_pre', default='_R', help="Paired-end prefix. String to insert between the file basename and the pair number when searching for read files. If your FASTQ files are names as [reads_R1.fastq, reads_R2.fastq] then use '_R1', or if reads_1.fastq then use '_1'. This option is only used when in paired end mode. default: _R1")
    parser.add_argument('--sample_postfix', default="", help="Postfix to append when looking for read files (ex lane number: '_001')")
    
    parser.add_argument('--trim', action='store_true', help="Use trimmomatic to perform adapter clipping.")
    parser.add_argument('--assembler', action='store', default='stringtie', choices=['stringtie', 'cufflinks'], help="Specify with transcript assembly program to use.")
    
    performance_group = add_runner_args(parser)
    #performance_group.add_argument('--skipalign', action='store_true', help="Skip the alignment process and only run the QC routines. Assumes you have previously aligned files in the proper locations.")
    ckpts = ['post_trim', 'pre_align', 'post_align', 'pre_assembly', 'post_assembly', 'pre_merge', 'post_merge']
    performance_group.add_argument('--resume', action='store', default=None, choices=ckpts, help='Resume the pipeline from this checkpoint.')

    args, additional_args = parser.parse_known_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices



    sys.stdout.write('Reading sample manifest.....\n')
    #sample manifest should be in the following TAB separated format (with headers):
    #Path    Basename    PE
    #/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true
    sample_manifest = pd.read_csv(args.manifest, sep='\t', comment='#', skip_blank_lines=True, true_values=['true', 'True', 'TRUE', '1'], false_values=['false', 'False', 'FALSE', '0'])
    samples = [Tuxedo2PipelineRawSample(s, args.genome, os.path.join(args.dest, s['Basename']), args.pe_pre, args.sample_postfix) for s in sample_manifest.to_dict(orient='records')]
    sys.stdout.write(' -> Found {count} item{plural} for processing.....\n'.format(count=len(samples), plural=('s' if len(samples) > 1 else '')))


    #get and run the read alignment pipeline
    pipeline = make_read_alignment_pipeline(args, additional_args)
    if args.resume is None or ckpts.index(args.resume) < ckpts.index('post_assembly'):
        pipeline.offset = args.resume
        runner = get_configured_runner(args, pipeline)
        runner.run(samples)
        sys.stdout.write("Completed alignment and initial quantification phase for all manifest items!\n")
        sys.stdout.write("=========================================================\n\n")
        sys.stdout.flush()
    
    
    #process samples from previous step, and generate new pseudo-sample for transcript merging
    merge_sample = Tuxedo2PipelineMergeSample({'Basename': 'TranscriptMergePseudoSample', 'Genome': args.genome, 'Dest': args.dest})
    for sample in samples:
        gtf = sample.find_files(lambda f: f.cxt.role == 'assembled_transcripts')
        sys.stdout.write(str(gtf))
        merge_sample.add_file(FileInfo(gtf[0].fullpath, FileContext.from_origin(sample.name)))
        
    pipeline = make_transcript_merge_pipeline(args)
    if args.resume is None or ckpts.index(args.resume) < ckpts.index('post_merge'):
        if ckpts.index(args.resume) >= ckpts.index('pre_merge'):
            pipeline.offset = args.resume
        runner = get_configured_runner(args, pipeline)
        runner.run([merge_sample])            
        sys.stdout.write("Completed merge of transcript assemblies from all samples!\n")
        sys.stdout.write("=========================================================\n\n")
        sys.stdout.flush()
    
    
    #re-process original samples, now using the merged transcripts
    merged_gtf = merge_sample.find_files(lambda f: f.cxt.role == 'merged_transcript_assembly')[0]
    for sample in samples:
        sample.add_file(FileInfo(merged_gtf, FileContext.from_origin('merged_transcript_assembly')))
        
    pipeline = make_transcript_quant_pipeline(args)
    runner = get_configured_runner(args, pipeline)
    runner.run(samples)
    sys.stdout.write("Completed re-quantification using merged transcript assembly for all manifest items!\n")
    sys.stdout.write("=========================================================\n\n")
    sys.stdout.flush()
    
#end main()

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
    
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('post_trim'))
    ################################
    
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
        
    
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('pre_align'))
    ################################
    
    
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
    
    def mapped_bam_resolver(cxt):
        return cxt.sample.find_files(lambda f: f.basename == '{}.bam'.format(cxt.sample.name))[0]
       
       
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('post_align'))
    ################################
    
    #QC functions:
    if (args.qc is not None) and (len(args.qc) > 0):
        
        if 'ism' in args.qc:
            from ThackTech.Pipelines.PipelineModules import InsertSizeMetrics
            x = InsertSizeMetrics.InsertSizeMetrics()
            x.set_resolver('bam', mapped_bam_resolver)
            pipeline.append_module(x)
    
    
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('pre_assembly'))
    ################################
            
    
    if args.assembler == 'cufflinks':
        # @todo: implement cufflinks support
        # 1. run cufflinks
        from ThackTech.Pipelines.PipelineModules import Cufflinks
        x = Cufflinks.Cufflinks(processors=args.threads)
        x.set_resolver('alignments', mapped_bam_resolver)
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.genome.genes_gtf)
        pipeline.append_module(x, critical=True)
        
        # 2. run cuffcompare? optional, I think....
        
    else:
        #for each sample, assemble transcripts
        #ex: stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 ERR188044_chrX.bam
        from ThackTech.Pipelines.PipelineModules import StringTie
        x = StringTie.StringTieQuant(processors=args.threads)
        x.set_resolver('alignments', mapped_bam_resolver)
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.genome.genes_gtf)  
        pipeline.append_module(x, critical=True)
        
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('post_assembly'))
    ################################

    return pipeline
#end make_read_alignment_pipeline()
    
   
    
    
def make_transcript_merge_pipeline(args):
    pipeline = AnalysisPipeline('Transcript Merge')
    
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('pre_merge'))
    ################################
    
    if args.assembler == 'cufflinks':
        # @todo: implement cufflinks support
        # 1. run cuffmerge
        from ThackTech.Pipelines.PipelineModules import Cuffmerge
        x = Cuffmerge.Cuffmerge(processors=args.threads)
        x.set_resolver('assemblies', lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.gtf'))
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.genome.genes_gtf)  
        pipeline.append_module(x, critical=True)
        
    else:
        # Merge Transcript assemblies for all samples
        from ThackTech.Pipelines.PipelineModules import StringTie
        x = StringTie.StringTieMerge(processors=args.threads)
        x.set_resolver('assemblies', lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.gtf'))
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.genome.genes_gtf)  
        pipeline.append_module(x, critical=True)
        
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('post_merge'))
    ################################

    
    #optional: Examine how the merged transcripts compare with the reference annotation
    #ex: gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf
    #from ThackTech.Pipelines.PipelineModules import GffCompare
    #x = GffCompare.GffCompare(processors=args.threads)
    #x.set_resolver('reference_gtf', lambda cxt: cxt.sample.genome.genes_gtf)
    #x.set_resolver('gtfs', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role == 'merged_transcript_assembly'))
    
    return pipeline
#end make_transcript_merge_pipeline()

   
    
def make_transcript_quant_pipeline(args):
    
    pipeline = AnalysisPipeline('Quantify Transcripts')
    
    def mapped_bam_resolver(cxt):
        return cxt.sample.find_files(lambda f: f.basename == '{}.bam'.format(cxt.sample.name))[0]
    
    if args.assembler == 'cufflinks':
        # @todo: implement cufflinks support
        # 1. run cuffquant
        from ThackTech.Pipelines.PipelineModules import Cuffquant
        x = Cuffquant.Cuffquant(processors=args.threads)
        x.set_resolver('alignments', mapped_bam_resolver)
        # @todo: should be merged gtf, think this is OK, but need to check!!!!!!!!
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role == 'merged_transcript_assembly'))  
        pipeline.append_module(x, critical=True)
        
        
        # 2. run cuffnorm? -> maybe need separate pipeline step.... needs all samples
        # 3. run cuffdiff? -> maybe need separate pipeline step.... needs all samples
        pass
        
    else:
        #Estimate transcript abundances and create table counts for Ballgown:
        #stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf ERR188044_chrX.bam
        from ThackTech.Pipelines.PipelineModules import StringTie
        x = StringTie.StringTieQuant(processors=args.threads)
        x.set_parameter('estimate', True)
        x.set_parameter('out_ballgown', True)
        x.set_resolver('alignments', mapped_bam_resolver)
        
        # @todo: should be merged gtf, think this is OK, but need to check!!!!!!!!
        x.set_resolver('guide_gff', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role == 'merged_transcript_assembly'))  
        pipeline.append_module(x, critical=True)
    
    return pipeline
#end make_transcript_quant_pipeline()







if __name__ == "__main__":
    main()










