#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, AnalysisPipelineCheckpoint, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner



#sample manifest should be in the following TAB separated format (with headers):
#Path    Basename    PE    Genome    Dest
#/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest

#BEFORE RUNNING SCRIPT WHEN USING LMOD
#module load java/1.8.0_121 samtools/0.1.19 intel/17.0.2 python/2.7.12 bedtools2/2.25.0 R-Project/3.3.3 bowtie2/2.2.9


class AlignmentPipelineSample(PipelineSample):

    def __init__(self, sample, pe_prefix='_R', postfix=""):
        super(AlignmentPipelineSample, self).__init__(sample['Basename'], sample['Genome'], sample['Dest'])
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
#end Pipelinesample


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', help="Manifest file containing sample information in tab separated format. Should contain the following columns (headers must be present): [Path], [Basename], [PE], [Genome], [Dest]")
    parser.add_argument('--bowtie-version', choices=['1', '2'], default='1', help="Version of bowtie to run")
    
    available_qc_choices = ['spp', 'pbc', 'ism', 'fpt', 'rpkm', 'fqscreen', 'fastqc']
    parser.add_argument('--qc', action='append', default=[], choices=available_qc_choices+['all'], help="Specify which QC pipelines to run after the alignment process completes. SPP is the cross-correlation analysis provided by ccQualityControl/phantompeakqualtools. PBC is the PRC bottleneck/library complexity estimation. ISM computes the distribution of insert size. FPT computes the \"BAM Fingerprint\" using DeepTools bamFingerprint program, and give a good idea of IP strength, especially for TF-like IPs. rpkm will generate a RPKM normalized BigWig from the aligned BAM file. All will run all available QC modules.")
    parser.add_argument('--pe_pre', default='_R', help="Paired-end prefix. String to insert between the file basename and the pair number when searching for read files. If your FASTQ files are names as [reads_R1.fastq, reads_R2.fastq] then use '_R1', or if reads_1.fastq then use '_1'. This option is only used when in paired end mode. default: _R1")
    parser.add_argument('--sample_postfix', default="", help="Postfix to append when looking for read files (ex lane number: '_001')")
    parser.add_argument('--unaligned', action='store_true', help='Output reads that fail to align to the reference genome.')
    parser.add_argument('--trim', action='store_true', help="Use trimmomatic to perform adapter clipping.")
    parser.add_argument('--override-dest', action='store', default=None, help="Override the destination read from the sample manifest.")
    
    performance_group = add_runner_args(parser)
    ckpts = ['post_trim', 'pre_align', 'post_align', 'post_pbc']
    performance_group.add_argument('--resume', action='store', default=None, choices=ckpts, help='Resume the pipeline from this checkpoint.')

    args, additional_args = parser.parse_known_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices

    sys.stdout.write('Reading sample manifest.....\n')
    #sample manifest should be in the following TAB separated format (with headers):
    #Path    Basename    PE    Genome    Dest
    #/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest
    sample_manifest = pd.read_csv(args.manifest, sep='\t', comment='#', skip_blank_lines=True, true_values=['true', 'True', 'TRUE', '1'], false_values=['false', 'False', 'FALSE', '0'])
    samples = [AlignmentPipelineSample(s, args.pe_pre, args.sample_postfix) for s in sample_manifest.to_dict(orient='records')]
    if args.override_dest is not None:
        sys.stdout.write("Override Destination is turned ON\n")
        sys.stdout.write('\t-> Setting destination for all samples to "{dest}"\n'.format(dest=os.path.abspath(args.override_dest)))
        for s in samples:
            s.dest = os.path.abspath(args.override_dest)
    sys.stdout.write('\t-> Found {count} item{plural} for processing.....\n'.format(count=len(samples), plural=('s' if len(samples) > 1 else '')))

    #get the pipeline
    pipeline = make_read_alignment_pipeline(args, additional_args)
    pipeline.offset = args.resume #handle any checkpoints
    
    #get the runner, then run the pipeline
    runner = get_configured_runner(args, pipeline)
    runner.run(samples) #runner blocks until processing is complete
    
    sys.stdout.write("Completed processing of all manifest items!\n")
    sys.stdout.write("=========================================================\n\n")
    sys.stdout.flush()
#end main()

def make_read_alignment_pipeline(args, additional_args):

    pipeline = AnalysisPipeline('Read Alignment and QC')
    
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferToShm
        x = TransferToShm.TransferToShm()
        x.set_parameter('shm_path', args.shm_path)
        pipeline.append_module(x, critical=True)
            
    #if not args.skipalign:
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
    
    
    if args.bowtie_version == '1':
        from ThackTech.Pipelines.PipelineModules import DecompressFiles
        x = DecompressFiles.DecompressFiles()
        x.set_resolver('files', resolve_bowtie1)
        pipeline.append_module(x, critical=True)
    
        from ThackTech.Pipelines.PipelineModules import BowtieAlign
        x = BowtieAlign.BowtieAlign()
        x.set_available_cpus(args.threads)
        x.set_parameter('unaligned', args.unaligned)
        x.set_parameter('additional_args', additional_args)
        x.set_resolver('fastq', resolve_bowtie1)
        pipeline.append_module(x, critical=True)
        
        from ThackTech.Pipelines.PipelineModules import RemoveDecompressedFiles
        x = RemoveDecompressedFiles.RemoveDecompressedFiles()
        pipeline.append_module(x, critical=True)
        
    elif args.bowtie_version == '2':#run bowtie2
        from ThackTech.Pipelines.PipelineModules import Bowtie2Align
        x = Bowtie2Align.Bowtie2Align()
        x.set_available_cpus(args.threads)
        x.set_parameter('unaligned', args.unaligned)
        x.set_parameter('additional_args', additional_args)
        x.set_resolver('fastq', resolve_bowtie1)
        pipeline.append_module(x, critical=True)
    else:
        raise ValueError("There is no bowtie version {}. Must be one of [1, 2].".format(args.bowtie_version))
        
    from ThackTech.Pipelines.PipelineModules import SamToBam
    x = SamToBam.SamToBam()
    x.set_available_cpus(args.threads)
    x.set_resolver('sam', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role == 'sam')[0])
    pipeline.append_module(x, critical=True)
    
    
    ################################
    # CHECKPOINT
    pipeline.append_module(AnalysisPipelineCheckpoint('post_align'))
    ################################
    
   
    if (args.qc is not None) and (len(args.qc) > 0):
        
        def qc_bt_bam_resolver(cxt):
            return cxt.sample.find_files(lambda f: f.basename == '{}.bam'.format(cxt.sample.name))[0]
        
        if 'pbc' in args.qc:
            from ThackTech.Pipelines.PipelineModules import PbcAnalysis
            x = PbcAnalysis.PbcAnalysis(processors=args.threads)
            x.set_resolver('bam', qc_bt_bam_resolver)
            pipeline.append_module(x)
            
        
        def qc_pbc_bam_resolver(cxt):
            return cxt.sample.find_files(lambda f: f.cxt.role == 'filtered_bam')[0]
        
        
        ################################
        # CHECKPOINT
        pipeline.append_module(AnalysisPipelineCheckpoint('post_pbc'))
        ################################
        
        #if 'spp' in args.qc:
        #    from ThackTech.Pipelines.PipelineModules import SPP
        #    pipeline.append_module(SPP.SPP())
            
        if 'ism' in args.qc:
            from ThackTech.Pipelines.PipelineModules import InsertSizeMetrics
            
            if 'pbc' in args.qc:
                x = InsertSizeMetrics.InsertSizeMetrics()
                x.set_resolver('bam', qc_pbc_bam_resolver)
                pipeline.append_module(x)
                
            x = InsertSizeMetrics.InsertSizeMetrics()
            x.set_resolver('bam', qc_bt_bam_resolver)
            pipeline.append_module(x)
            
        
        if 'fpt' in args.qc:
            from ThackTech.Pipelines.PipelineModules import BamFingerprint
            x = BamFingerprint.BamFingerprint(processors=args.threads)
            if 'pbc' in args.qc:
                x.set_resolver('bams', lambda cxt: [qc_bt_bam_resolver(cxt), qc_pbc_bam_resolver(cxt)])
            else:
                x.set_resolver('bams', lambda cxt: [qc_bt_bam_resolver(cxt)])
            pipeline.append_module(x)
            
            
        if 'rpkm' in args.qc:
            from ThackTech.Pipelines.PipelineModules import BamCoverage

            if 'pbc' in args.qc:
                x = BamCoverage.BamCoverage(processors=args.threads)
                x.set_resolver('bam', qc_pbc_bam_resolver)
                pipeline.append_module(x)
                
            x = BamCoverage.BamCoverage(processors=args.threads)
            x.set_resolver('bam', qc_bt_bam_resolver)
            pipeline.append_module(x)
            
            
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
    
    return pipeline
#end make_read_alignment_pipeline()



if __name__ == "__main__":
    main()









