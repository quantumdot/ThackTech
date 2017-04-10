#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner



#sample manifest should be in the following TAB separated format (with headers):
#Path    Basename    PE    Genome    Dest
#/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest

#BEFORE RUNNING SCRIPT WHEN USING LMOD
#module load java/1.8.0_121 samtools/0.1.19 intel/17.0.2 python/2.7.12 bedtools2/2.25.0 R-Project/3.3.3 bowtie2/2.2.9

gopts = {
    'spp_path': '/home/josh/scripts/phantompeakqualtools/run_spp.R',
    'shm_dir': '/mnt/ramdisk/btalign',
    'pe_prefix': '_R',
}

    
class AlignmentPipelineSample(PipelineSample):

    def __init__(self, sample):
        super(AlignmentPipelineSample, self).__init__(sample['Basename'], sample['Genome'], sample['Dest'])
        self.set_attribute('PE', ('PE' in sample and sample['PE']))
        self.discover_files(sample['Path'])
    #end __init__()
    
    def discover_files(self, path):
        files = []
        compressed_extensions = ['.gz', '.bz2', '.zip', '.tar', '.tar.gz']
        if self.get_attribute('PE'):
            base = os.path.join(path, self.name + gopts['pe_prefix'])
            #print "trying %s, %s" % (base+'1.fastq', base+'2.fastq')
            if os.path.exists(base+'1.fastq') and os.path.exists(base+'2.fastq'):
                files.append(FileInfo(base+'1.fastq', FileContext.from_origin('reads'), mate=1))
                files.append(FileInfo(base+'2.fastq', FileContext.from_origin('reads'), mate=2))
            else:
                for ext in compressed_extensions:
                    if os.path.exists(base+'1.fastq'+ext) and os.path.exists(base+'2.fastq'+ext):
                        files.append(FileInfo(base+'1.fastq'+ext, FileContext.from_origin('reads'), mate=1))
                        files.append(FileInfo(base+'2.fastq'+ext, FileContext.from_origin('reads'), mate=2))
                        continue
        else:
            base = os.path.join(path, self.name+'.fastq')
            if os.path.exists(base):
                files.append(FileInfo(base, FileContext.from_origin('reads')))
            else:
                for ext in compressed_extensions:
                    if os.path.exists(base+ext):
                        files.append(FileInfo(base+ext, FileContext.from_origin('reads')))
                        continue
        if len(files) <= 0:
            raise IOError('Unable to find '+('PE ' if self.get_attribute('PE') else '')+'reads for '+self.name)
        for f in files:
            self.add_file(f)
    #end find_files()
#end Pipelinesample


def main():
    #total_start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', help="Manifest file containing sample information in tab separated format. Should contain the following columns (headers must be present): [Path], [Basename], [PE], [Genome], [Dest]")
    parser.add_argument('--bowtie-version', choices=['1', '2'], default='1', help="Version of bowtie to run")
    
    available_qc_choices = ['spp', 'pbc', 'ism', 'fpt', 'rpkm', 'fqscreen']
    parser.add_argument('--qc', action='append', default=[], choices=available_qc_choices+['all'], help="Specify which QC pipelines to run after the alignment process completes. SPP is the cross-correlation analysis provided by ccQualityControl/phantompeakqualtools. PBC is the PRC bottleneck/library complexity estimation. ISM computes the distribution of insert size. FPT computes the \"BAM Fingerprint\" using DeepTools bamFingerprint program, and give a good idea of IP strength, especially for TF-like IPs. rpkm will generate a RPKM normalized BigWig from the aligned BAM file. All will run all available QC modules.")
    parser.add_argument('--pe_pre', default='_R', help="Paired-end prefix. String to insert between the file basename and the pair number when searching for read files. If your FASTQ files are names as [reads_R1.fastq, reads_R2.fastq] then use '_R1', or if reads_1.fastq then use '_1'. This option is only used when in paired end mode. default: _R1")
    parser.add_argument('--unaligned', action='store_true', help='Output reads that fail to align to the reference genome.')
    parser.add_argument('--trim', action='store_true', help="Use trimmomatic to perform adapter clipping.")
    
    performance_group = add_runner_args(parser)
    performance_group.add_argument('--skipalign', action='store_true', help="Skip the alignment process and only run the QC routines. Assumes you have previously aligned files in the proper locations.")
    
    
    args, additional_args = parser.parse_known_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices

    sys.stdout.write('Reading sample manifest.....\n')
    sample_manifest = pd.read_csv(args.manifest, sep='\t', comment='#', skip_blank_lines=True, true_values=['true', 'True', 'TRUE', '1'], false_values=['false', 'False', 'FALSE', '0'])
    samples = [AlignmentPipelineSample(s) for s in sample_manifest.to_dict(orient='records')]
    sys.stdout.write('\t-> Found %d item%s for processing.....\n' % (len(sample_manifest.index), ('s' if len(sample_manifest.index) > 1 else '')))

    #sample manifest should be in the following TAB separated format (with headers):
    #Path    Basename    PE    Genome    Dest
    #/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest

    gopts['pe_prefix'] = args.pe_pre

    
    pipeline = make_read_alignment_pipeline(args, additional_args)
    runner = get_configured_runner(args, pipeline)
    runner.run(samples)
    
    sys.stdout.write("Completed processing of all manifest items!\n")
    #sys.stdout.write("\t-> Processing entire manifest took %s\n\n" % (Common.human_time_diff(total_start_time, time.time()),))
    sys.stdout.write("=========================================================\n\n")
    sys.stdout.flush()
#end main()

def make_read_alignment_pipeline(args, additional_args):
    pipeline = AnalysisPipeline('Read Alignment')
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferToShm
        x = TransferToShm.TransferToShm()
        x.set_parameter('shm_path', args.shm_path)
        pipeline.append_module(x, critical=True)
            
    if not args.skipalign:
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
        
        if 'fqscreen' in args.qc:    
            from ThackTech.Pipelines.PipelineModules import FastqScreen
            x = FastqScreen.FastqScreen()
            x.set_resolver('fastqs', resolve_bowtie1)
            pipeline.append_module(x)
    
    
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
    
    
    else: #we are skipping alignment process
        from ThackTech.Pipelines.PipelineModules import ReadOutputManifest
        pipeline.append_module(ReadOutputManifest.ReadOutputManifest(), critical=True)
        
            
    if (args.qc is not None) and (len(args.qc) > 0):
        if 'pbc' in args.qc:
            from ThackTech.Pipelines.PipelineModules import PbcAnalysis
            x = PbcAnalysis.PbcAnalysis(processors=args.threads)
            x.set_resolver('bam', lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.bam')[0])
            pipeline.append_module(x)
            
            def qc_bam_resolver(cxt):
                return cxt.sample.find_files(lambda f: f.cxt.role == 'filtered_deduplicated_bam')[0]
        else:
            def qc_bam_resolver(cxt):
                return cxt.sample.find_files(lambda f: f.ext == '.bam')[0]
        
        #if 'spp' in args.qc:
        #    from ThackTech.Pipelines.PipelineModules import SPP
        #    pipeline.append_module(SPP.SPP())
            
        if 'ism' in args.qc:
            from ThackTech.Pipelines.PipelineModules import InsertSizeMetrics
            x = InsertSizeMetrics.InsertSizeMetrics()
            x.set_resolver('bam', qc_bam_resolver)
            pipeline.append_module(x)
        
        if 'fpt' in args.qc:
            from ThackTech.Pipelines.PipelineModules import BamFingerprint
            x = BamFingerprint.BamFingerprint(processors=args.threads)
            x.set_resolver('bams', lambda cxt: cxt.sample.find_files(lambda f: f.basename in ['{}.bam'.format(cxt.sample.name), '{}.filt.srt.nodup.bam'.format(cxt.sample.name)]))
            pipeline.append_module(x)
            
        if 'rpkm' in args.qc:
            from ThackTech.Pipelines.PipelineModules import BamToRpkmNormBigWig
            x = BamToRpkmNormBigWig.BamToRpkmNormBigWig(processors=args.threads)
            x.set_resolver('bam', qc_bam_resolver)
            pipeline.append_module(x)
    
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
    
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    pipeline.append_module(OutputManifest.OutputManifest())
    
    return pipeline
#end make_read_alignment_pipeline()



if __name__ == "__main__":
    main()









