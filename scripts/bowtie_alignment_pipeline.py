#!/usr/bin/env python

import os
import multiprocessing
import pandas as pd
import argparse
import sys
from ThackTech import Common
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline


#sample manifest should be in the following TAB separated format (with headers):
#Path    Basename    PE    Genome    Dest
#/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest

gopts = {
    'spp_path': '/home/josh/scripts/phantompeakqualtools/run_spp.R',
    'shm_dir': '/mnt/ramdisk/btalign',
    'pe_prefix': '_R',
}
try:
    cpu_count = multiprocessing.cpu_count()
except:
    cpu_count = 1
if cpu_count is None:
    cpu_count = 1

    
class AlignmentPipelineSample(PipelineSample):

    def __init__(self, sample):
        PipelineSample.__init__(self, sample['Basename'], sample['Genome'], sample['Dest'])
        self.add_attribute('PE', ('PE' in sample and sample['PE']))
        self.discover_files(sample['Path'])
    #end __init__()
    
    def discover_files(self, path):
        files = []
        is_compressed = False
        compressed_extensions = ['.gz', '.bz2', '.zip', '.tar', '.tar.gz']
        if self.get_attribute('PE'):
            base = os.path.join(path, self.name + gopts['pe_prefix'])
            #print "trying %s, %s" % (base+'1.fastq', base+'2.fastq')
            if os.path.exists(base+'1.fastq') and os.path.exists(base+'2.fastq'):
                files.append(base+'1.fastq')
                files.append(base+'2.fastq')
            else:
                for ext in compressed_extensions:
                    if os.path.exists(base+'1.fastq'+ext) and os.path.exists(base+'2.fastq'+ext):
                        files.append(base+'1.fastq'+ext)
                        files.append(base+'2.fastq'+ext)
                        continue
                        is_compressed = True
        else:
            base = os.path.join(path, self.name+'.fastq')
            if os.path.exists(base):
                files.append(base)
            else:
                for ext in compressed_extensions:
                    if os.path.exists(base+ext):
                        files.append(base+ext)
                        continue
                        is_compressed = True
        if len(files) <= 0:
            raise IOError('Unable to find '+('PE ' if self.get_attribute('PE') else '')+'reads for '+self.name)
        self.add_file('source', 'fastq', files)
    #end find_files()
#end Pipelinesample


def main():
    #total_start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', help="Manifest file containing sample information in tab separated format. Should contain the following columns (headers must be present): [Path], [Basename], [PE], [Genome], [Dest]")
    parser.add_argument('--bowtie-version', choices=['1', '2'], default='1', help="Version of bowtie to run")
    
    available_qc_choices = ['spp', 'pbc', 'ism', 'fpt', 'rpkm']
    parser.add_argument('--qc', action='append', default=[], choices=available_qc_choices+['all'], help="Specify which QC pipelines to run after the alignment process completes. SPP is the cross-correlation analysis provided by ccQualityControl/phantompeakqualtools. PBC is the PRC bottleneck/library complexity estimation. ISM computes the distribution of insert size. FPT computes the \"BAM Fingerprint\" using DeepTools bamFingerprint program, and give a good idea of IP strength, especially for TF-like IPs. rpkm will generate a RPKM normalized BigWig from the aligned BAM file. All will run all available QC modules.")
    parser.add_argument('--pe_pre', default='_R', help="Paired-end prefix. String to insert between the file basename and the pair number when searching for read files. If your FASTQ files are names as [reads_R1.fastq, reads_R2.fastq] then use '_R1', or if reads_1.fastq then use '_1'. This option is only used when in paired end mode. default: _R1")
    parser.add_argument('--unaligned', action='store_true', help='Output reads that fail to align to the reference genome.')
    parser.add_argument('--trim', action='store_true', help="Use trimmomatic to perform adapter clipping.")
    
    performance_group = parser.add_argument_group('Performance')
    performance_group.add_argument('-p', '--threads', type=int, default=cpu_count, help="Number of processors to use for processing.")
    #performance_group.add_argument('--noparallel', action='store_true', help='Completly turn off parallelization of tasks (no use of the thread pool).')
    performance_group.add_argument('--shm', action='store_true', help="Use ramfs for file IO.")
    performance_group.add_argument('--shm-path', action='store', default='/mnt/ramdisk/bowtie'+'_'+str( os.getuid() ), help='When --shm is passed, the path to use for ram disk storage. Individual samples will have dedicated subfolders on this path. Please ensure this path has appropriate permissions.')
    performance_group.add_argument('--skipalign', action='store_true', help="Skip the alignment process and only run the QC routines. Assumes you have previously aligned files in the proper locations.")
    performance_group.add_argument('--runner', action='store', default='parallel', choices=['slurm', 'parallel', 'serial'], help="Which pipeline runner to use.")
    
    args, additional_args = parser.parse_known_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices

    sys.stdout.write('Reading sample manifest.....\n')
    sample_manifest = pd.read_csv(args.manifest, sep='\t', comment='#', skip_blank_lines=True, true_values=['true', 'True', 'TRUE', '1'], false_values=['false', 'False', 'FALSE', '0'])
    samples = [AlignmentPipelinesample(s) for s in sample_manifest.to_dict(orient='records')]
    sample_count = len(samples)
    sys.stdout.write('\t-> Found %d item%s for processing.....\n' % (len(sample_manifest.index), ('s' if len(sample_manifest.index) > 1 else '')))

    #sample manifest should be in the following TAB separated format (with headers):
    #Path    Basename    PE    Genome    Dest
    #/path/to/fastq    anti_H3K18Ac_K562_WCE_CAGATC_ALL    true    hg19    /path/to/bam/dest

    gopts['pe_prefix'] = args.pe_pre

    
    pipeline = make_read_alignment_pipeline(args, additional_args)
    #sys.stdout.write(pipeline.documentation())
    
    
    #run the pipeline!
    if args.runner == 'slurm':
        from ThackTech.Pipelines import SlurmPipelineRunner
        runner = SlurmPipelineRunner(pipeline, partition="main", nodes=1, threads_per_node=args.threads, time_limit="5:00:00")
    elif args.runner == 'parallel':
        from ThackTech.Pipelines import ParallelPipelineRunner
        runner = ParallelPipelineRunner(pipeline, args.threads)
    else: #serial runner
        from ThackTech.Pipelines import SerialPipelineRunner
        runner = SerialPipelineRunner(pipeline)
    
    
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
            x = Trimmomatic.Trimmomatic()
            x.set_available_cpus(args.threads)
            pipeline.append_module(x, critical=True)
    
    
        if args.bowtie_version == '1':
            from ThackTech.Pipelines.PipelineModules import DecompressFiles
            x = DecompressFiles.DecompressFiles()
            x.set_resolver('files', lambda s: [['source', 'fastq']])
            pipeline.append_module(x, critical=True)
        
            from ThackTech.Pipelines.PipelineModules import BowtieAlign
            x = BowtieAlign.BowtieAlign()
            x.set_available_cpus(args.threads)
            x.set_parameter('unaligned', args.unaligned)
            x.set_resolver('fastq', lambda s: s.get_file('source', 'fastq'))
            x.set_parameter('additional_args', additional_args)
            pipeline.append_module(x, critical=True)
            
            from ThackTech.Pipelines.PipelineModules import SamToBam
            x = SamToBam.SamToBam()
            x.set_available_cpus(args.threads)
            x.set_resolver('sam', lambda s: s.get_file('BowtieAlign', 'sam'))
            pipeline.append_module(x, critical=True)
            
            from ThackTech.Pipelines.PipelineModules import RemoveDecompressedFiles
            x = RemoveDecompressedFiles.RemoveDecompressedFiles()
            pipeline.append_module(x, critical=True)
            
        else:#run bowtie2
            from ThackTech.Pipelines.PipelineModules import Bowtie2Align
            x = Bowtie2Align.Bowtie2Align()
            x.set_available_cpus(args.threads)
            #x.set_parameter('duplicates', args.duplicates)
            #x.set_parameter('bw', args.bw)
            pipeline.append_module(x, critical=True)
            
            from ThackTech.Pipelines.PipelineModules import SamToBam
            x = SamToBam.SamToBam()
            x.set_available_cpus(args.threads)
            x.set_resolver('sam', lambda s: s.get_file('Bowtie2Align', 'sam'))
            pipeline.append_module(x, critical=True)
    
    
    else: #we are skipping alignment process
        from ThackTech.Pipelines.PipelineModules import ReadOutputManifest
        pipeline.append_module(ReadOutputManifest.ReadOutputManifest(), critical=True)
        
            
    if (args.qc is not None) and (len(args.qc) > 0):
        if 'pbc' in args.qc:
            from ThackTech.Pipelines.PipelineModules import PbcAnalysis
            x = PbcAnalysis.PbcAnalysis()
            x.set_resolver('bam', lambda s: s.get_file())
            pipeline.append_module()
        
        if 'spp' in args.qc:
            from ThackTech.Pipelines.PipelineModules import SPP
            pipeline.append_module(SPP.SPP())
            
        if 'ism' in args.qc:
            from ThackTech.Pipelines.PipelineModules import InsertSizeMetrics
            pipeline.append_module(InsertSizeMetrics.InsertSizeMetrics())
        
        if 'fpt' in args.qc:
            from ThackTech.Pipelines.PipelineModules import BamFingerprint
            pipeline.append_module(BamFingerprint.BamFingerprint())
            
        if 'rpkm' in args.qc:
            from ThackTech.Pipelines.PipelineModules import RPKMNormBigWig
            x = RPKMNormBigWig.RPKMNormBigWig()
            x.set_available_cpus(args.threads)
            x.set_resolver('bam', lambda s: s.get_file('SamToBam', 'bam'))
            pipeline.append_module(x)
    
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
    
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    pipeline.append_module(OutputManifest.OutputManifest())
    
    return pipeline
#end make_read_alignment_pipeline()
    
    
    
    
    
    
    
    
    
    
    # Common.ensure_dir(gopts['shm_dir'])
        
        
        
    
    # if not args.skipalign:
        # show_bowtie_version()
        # sys.stdout.write('\n=========================================================\n\n')
        # for sidx, sample in sample_manifest.iterrows():
            # sample_start_time = time.time()
            # sys.stdout.write("Processing sample %s......\n" % (sample['Basename'],))
            # sys.stdout.write("\t-> Wall clock: %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"),))
            # sys.stdout.flush()
            
            # working_dest = gopts['shm_dir'] if args.shm else sample['Dest']
            # Common.ensure_dir(sample['Dest'])

            # #find the read files
            # read_files, iscompressed = find_files(sample['Path'], sample['Basename'], sample['PE'])
            
            # #perform trimming
            # if args.trim:
                # read_files = run_trimmomatic(sample['Basename'], read_files, sample['PE'], os.path.join(sample['Dest'], 'trimmedfq'), args.threads)
                # iscompressed = True #trimmomatic is setup here to output gzip files
            
            # #handle file compression and use of shm
            # if iscompressed:
                # sys.stdout.write("\t-> Decompressing reads....\n")
                # sys.stdout.flush()
                # #if we are using shm, decompress to shm, else decompress to the src directory
                # read_files = decompress_reads(read_files, gopts['shm_dir'] if args.shm else sample['Path'])
            # elif args.shm:
                # sys.stdout.write("\t-> Moving reads to RAMFS....\n")
                # sys.stdout.flush()
                # #the reads are already uncompressed, so lets just copy them to shm
                # read_files = cp_to_shm(read_files)


            

            # #lets do a bit of cleanup
            # if iscompressed or args.shm:
                # #remove files that we decompressed or copied to shm, as they are no longer needed
                # for f in read_files:
                    # os.remove(f)

            # #Now lets convert SAM to BAM
            # bam_dest = sam_to_bam(sam_destination, args.threads)
            
            # #if we were working in shm, move the final final output to the dest directory
            # if args.shm:
                # sys.stdout.write("\t-> Moving output from RAMFS to final destination...\n")
                # sys.stdout.flush()
                # files_to_move = [bam_dest, bam_dest+'.bai']
                # if args.unaligned:
                    # unaln_files = []
                    # if sample['PE']:
                        # unaln_files.append(os.path.join(working_dest, sample['Basename']+'_unaligned_1.fastq'))
                        # unaln_files.append(os.path.join(working_dest, sample['Basename']+'_unaligned_2.fastq'))
                    # else:
                        # unaln_files.append(os.path.join(working_dest, sample['Basename']+'_unaligned.fastq'))
                    # proc = subprocess.Popen(['gzip'] + unaln_files)
                    # proc.communicate()
                    # for f in unaln_files:
                        # files_to_move.append(f+'.gz')
                # proc = subprocess.Popen(['mv'] + files_to_move + [sample['Dest']])
                # proc.communicate()
                # #shutil.move(bam_dest, os.path.join(sample['Dest'], os.path.basename(bam_dest)))
                # #shutil.move(bam_dest+'.bai', os.path.join(sample['Dest'], os.path.basename(bam_dest)+'.bai')) #dont forget about the index!!!

            # #let the user know our progress
            # sys.stdout.write("Done processing sample %s!\n" % (sample['Basename'],))
            # sys.stdout.write('\t-> See output at \"%s\"\n' % (sample['Dest'],))
            # sys.stdout.write("\t-> Processing sample took %s\n" % (Common.human_time_diff(sample_start_time, time.time()),))
            # sys.stdout.write("=========================================================\n\n")
            # sys.stdout.flush()
    # #end if not args.skipalign
    
    
    # if (args.qc is not None) and (len(args.qc) > 0):
        # sys.stdout.write("Preparing to run QC assays.....\n")
    
        # for sidx, sample in sample_manifest.iterrows():
            # sys.stdout.write('Processing sample "%s"....\n' % (sample['Basename'],))
            # sys.stdout.flush()
            # bam = os.path.join(sample['Dest'], sample['Basename']+'.bam')
            # if 'pbc' in args.qc:
                
            
            # if 'spp' in args.qc:
                
                
            # if 'ism' in args.qc:
                
            
            # if 'fpt' in args.qc:
                # sys.stdout.write('\t-> Running BAM fingerprint analysis...\n')
                # sys.stdout.flush()
                # fpt_dir = os.path.join(sample['Dest'], 'fingerprint')
                # Common.ensure_dir(fpt_dir)
                # bam_fingerprint_args = [
                    # 'bamFingerprint',
                    # '--bamfiles', bam,
                    # '--labels', sample['Basename'],
                    # '--plotFile', os.path.join(fpt_dir, sample['Basename']+'.fingerprint.pdf'),
                    # '--plotFileFormat', 'pdf',
                    # '--outRawCounts', os.path.join(fpt_dir, sample['Basename']+'.fingerprint_counts.txt'),
                    # '--numberOfProcessors', str(args.threads),
                    # '--binSize', '500', #default
                    # '--numberOfsamples', '1000000', #2x default (default is 0.5e6)
                    # #'--verbose'
                # ]
                # sys.stdout.write("\t-> BAM fingerprint analysis......")
                # sys.stdout.write("\n..............................................\n")
                # sys.stdout.write(" ".join(bam_fingerprint_args))
                # sys.stdout.write("\n..............................................\n")
                # sys.stdout.flush()
                # proc = subprocess.Popen(bam_fingerprint_args)
                # proc.communicate()
                # sys.stdout.write('\t-> Completed fingerprint analysis...\n')
                # sys.stdout.flush()
                
                
            # sys.stdout.write('Completed QC assays for sample "%s"....\n' % (sample['Basename'],))

    # sys.stdout.write("Completed processing of all manifest items!\n")
    # sys.stdout.write("\t-> Processing entire manifest took %s\n\n" % (Common.human_time_diff(total_start_time, time.time()),))
    # sys.stdout.write("=========================================================\n\n")
    # sys.stdout.flush()
#end main()


if __name__ == "__main__":
    main()









