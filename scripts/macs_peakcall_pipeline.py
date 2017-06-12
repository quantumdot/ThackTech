import os
import pandas as pd
import argparse
import sys
import copy
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline, FileInfo, FileContext
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner
from ThackTech.Pipelines.FileInfo import FileInfo, FileContext

#BEFORE RUNNING SCRIPT WHEN USING LMOD
#module load java/1.8.0_121 samtools/0.1.19 intel/17.0.2 python/2.7.12 bedtools2/2.25.0 R-Project/3.3.3 bowtie2/2.2.9


#manifest should be structured as follows (with headers)
#Name    Group    Treatment    Control    Genome    Broad    Dest


class MacsPipelineSample(PipelineSample):

    def __init__(self, sample):
        PipelineSample.__init__(self, sample['Name'], sample['Genome'], sample['Dest'])
        self.add_file(FileInfo(sample['Treatment'], FileContext.from_origin('treatment')))
                      
        if 'Control' in sample:
            self.add_file(FileInfo(sample['Control'], FileContext.from_origin('control')))
        if 'Broad' in sample and sample['Broad']:
            self.set_attribute('broad', True)
        if 'Group' in sample:
            self.set_attribute('group', sample['Group'])
    #end __init__()
#end PipelineSample




def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Given a sample manifest, runs MACS.',
        #epilog='The installed version of MACS is: '+show_macs_version(gopts['macs1_path'], None, False)+'; '+show_macs_version(gopts['macs2_path'], None, False)
    )
    parser.add_argument('manifest', help="Manifest file containing sample information. Manifest file should contain the following columns (tsv, case sensitive): [Name] [Genome] [Dest] [Treatment] [Control<optional>] [Broad<optional; bool-like>]")
    parser.add_argument('--macs-version', choices=['macs1', 'macs2'], help="Specifies which version of MACS to run")
    parser.add_argument('-d', '--duplicates', default='auto', help="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}.")
    #parser.add_argument('-f', '--format', default=None, choices=["BED", "SAM", "BAM", "BAMPE", "BOWTIE", "ELAND", "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET"], help="Format of the alignment files for the entire run. If not specified, format will be auto-detected.")
    parser.add_argument('-bw', default=300, type=int, help="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab.")
    parser.add_argument('--sigout', default='bdg', choices=['wig', 'bdg'], help="Output type for signal. Either wig or bedgraph.")
    available_qc_choices = ['chance', 'fingerprint', 'frip', 'rpkmbw']
    parser.add_argument('--qc', choices=available_qc_choices+['all'], default=[], action='append', help="Specify which QC routines to run. chance runs the IPStrength and spectrum modules from the chance package (Song et al. Genome Biology 2012. doi:10.1186/gb-2012-13-10-r98). fingerprint runs the DeepTools BamFingerprint module. FRiP computes the Fraction of Reads in Peaks. all will run all available QC modules. rpkmbw generates a RPKM normalized BigWig.")
    parser.add_argument('--idr', action='store_true', help="Perform IDR analysis. Requires that the manifest contain an additional [Group] column that describes group names for replicates. Will generate pooled and psuedoreplicates and pass them through the same processing pipeline as the explicitly defined manifest items. IDR is only supported in combination with MACS2.")
    parser.add_argument('--ignore-control', action='store_true', help="Ignore control data present in the sample manifest.")
    
    performance_group = add_runner_args(parser)
    performance_group.add_argument('--skipmacs', action='store_true', help="Skip the peak calling process and only run the QC routines. Assumes files are in the location specified by the manifest.")
    
    args = parser.parse_args()
    
    if 'all' in args.qc:
        args.qc = available_qc_choices
        
    if args.idr and args.macs_version == 'macs1':
        sys.stdout.write("IDR analysis is only supported in combination with MACS2!\n")
        sys.stdout.write("Exiting.....\n\n")
        sys.stdout.flush()
        return

    sys.stdout.write('Reading sample manifest.....\n')
    sample_manifest = pd.read_csv(args.manifest, sep='\t', comment='#', skip_blank_lines=True, true_values=['true', 'True', 'TRUE', '1'], false_values=['false', 'False', 'FALSE', '0'])
    
    if args.ignore_control:
        sample_manifest.drop('Control', axis=1, inplace=True)
    
    samples = [MacsPipelineSample(s) for s in sample_manifest.to_dict(orient='records')]
    sample_count = len(samples)
    sys.stdout.write('-> Found %d item%s for processing.\n\n' % (sample_count, ('s' if sample_count > 1 else '')))
    
    
    if args.idr:
        #first generate the pooled sample.
        generate_replicate_pools(samples, args)
        
        #generate pseudoreplicates
        generate_psuedoreplicate_samples(samples, args)
    
    ##
    # Perform Peak calling on samples
    ##
    run_peakcalling(samples, args)
    
    
    if args.idr:
        #run the IDR analysis
        run_idr_analysis(samples, args)
#end main()

def run_pipeline(pipeline, samples, args):
    runner = get_configured_runner(args, pipeline)
    runner.run(samples)
#end get_runner()

def run_read_manifest_pipeline(samples):
    pipeline = AnalysisPipeline('Reading Sample Manifest')
    from ThackTech.Pipelines.PipelineModules import ReadOutputManifest
    pipeline.append_module(ReadOutputManifest.ReadOutputManifest(), critical=True)
    
    from ThackTech.Pipelines import SerialPipelineRunner
    runner = SerialPipelineRunner(pipeline)
    runner.run(samples)
#end run_read_manifest_pipeline()

def run_peakcalling(samples, args):
    main_pipeline_samples_to_run = []
    main_pipeline_samples_to_skip = []
    for s in samples:
        out_manifest_dest = os.path.join(s.dest, s.name+'_output_manifest.tsv')
        if os.path.isfile(out_manifest_dest) and not args.skipmacs:
            sys.stdout.write('-> It appears that %s already has been completed. Skipping this sample...\n' % (s.name,))
            main_pipeline_samples_to_skip.append(s)
        else:
            main_pipeline_samples_to_run.append(s)
    sys.stdout.flush()
    
    run_read_manifest_pipeline(main_pipeline_samples_to_skip)
    run_pipeline(make_peak_calling_and_qc_pipeline(args), main_pipeline_samples_to_run, args)
            
    sys.stdout.write("Completed processing of all manifest items!\n\n")
    sys.stdout.write("--------------------------------------------\n\n")
    samples = (main_pipeline_samples_to_run + main_pipeline_samples_to_skip)
#end run_peakcalling()

def make_peak_calling_and_qc_pipeline(args):
    pipeline = AnalysisPipeline('Peak Calling and QC')
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferToShm
        x = TransferToShm.TransferToShm()
        x.set_parameter('shm_path', args.shm_path)
        pipeline.append_module(x, critical=True)
        
    if not args.skipmacs: 
        if args.macs_version == 'macs2':
            from ThackTech.Pipelines.PipelineModules import MACS2Peakcall
            x = MACS2Peakcall.MACS2Peakcall(processors=args.threads)
            x.set_parameter('duplicates', args.duplicates)
            x.set_parameter('bandwith', args.bw)
            x.set_resolver('treatments', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.cxt.role == 'treatment'))
            x.set_resolver('controls', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.cxt.role == 'control'))
            pipeline.append_module(x, critical=True)
        else:
            from ThackTech.Pipelines.PipelineModules import MACS1Peakcall
            x = MACS1Peakcall.MACS1Peakcall(processors=args.threads)
            x.set_parameter('duplicates', args.duplicates)
            x.set_parameter('bandwith', args.bw)
            x.set_parameter('sigout', args.sigout)
            x.set_resolver('treatments', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.cxt.role == 'treatment'))
            x.set_resolver('controls', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.cxt.role == 'control'))
            pipeline.append_module(x, critical=True)
            
        
    from ThackTech.Pipelines.PipelineModules import GeneratePlainBed
    x = GeneratePlainBed.GeneratePlainBed(processors=args.threads)
    x.set_resolver('encode_beds', lambda cxt: cxt.sample.find_files(lambda f: f.ext in ['.narrowPeak', '.broadPeak', '.gappedPeak']))
    pipeline.append_module(x)
    
    from ThackTech.Pipelines.PipelineModules import ConvertBedgraphToBigWig
    x = ConvertBedgraphToBigWig.ConvertBedgraphToBigWig(processors=args.threads)
    x.set_resolver('bedgraphs', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.role in ['treatment_signal', 'control_signal'] and f.ext == '.bdg'))
    pipeline.append_module(x)
        
    #if 'chance' in args.qc:
    #    from ThackTech.Pipelines.PipelineModules import ChanceAnalysis
    #    pipeline.append_module(ChanceAnalysis.ChanceAnalysis())
        
    if 'fingerprint' in args.qc:
        from ThackTech.Pipelines.PipelineModules import BamFingerprint
        x = BamFingerprint.BamFingerprint(processors=args.threads)
        x.set_resolver('bams', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.ext == '.bam'))
        pipeline.append_module(x)
        
    if 'frip' in args.qc:
        from ThackTech.Pipelines.PipelineModules import FRiPAnalysis
        x = FRiPAnalysis.FRiPAnalysis(processors=args.threads)
        def resolve_frip_input(cxt):
            finders = [
                lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.narrowPeak'),
                lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.broadPeak'),
                lambda cxt: cxt.sample.find_files(lambda f: f.ext == '.bed')
            ]
            for finder in finders:
                results = finder(cxt)
                if results is not None and len(results) > 0:
                    for r in results:
                        if 'peaks' in r.cxt.role.lower():
                            return r
            return None
        #end resolve_frip_input()
        x.set_resolver('bed', resolve_frip_input)
        x.set_resolver('bams', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.ext == '.bam'))
        pipeline.append_module(x)
        
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
        
    if 'rpkmbw' in args.qc:
        from ThackTech.Pipelines.PipelineModules import BamToRpkmNormBigWig
        x = BamToRpkmNormBigWig.BamToRpkmNormBigWig(name="Treatment_RPKM_Norm", processors=args.threads)
        x.set_resolver('bam',  lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.ext == '.bam' and f.cxt.role == 'treatment'))
        pipeline.append_module(x) 
        
        x = BamToRpkmNormBigWig.BamToRpkmNormBigWig(name="Control_RPKM_Norm", processors=args.threads)
        x.set_resolver('bam', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin and f.ext == '.bam' and f.cxt.role == 'control'))
        pipeline.append_module(x)    
        
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    pipeline.append_module(OutputManifest.OutputManifest())
    
    return pipeline
#end make_peak_calling_and_qc_pipeline


def generate_replicate_pools(samples, args):
    group_samples = {}
    for sample in samples:
        group = sample.get_attribute('group')
        if group not in group_samples:
            group_samples[group] = PipelineSample(group, sample.genome, sample.dest)
            group_samples[group].set_attribute('template', sample.name)
            
        for sf in sample.find_files(lambda f: f.cxt.is_origin()):
            c = len(group_samples[group].find_files(lambda f: f.cxt.is_origin() and f.cxt.role == sf.role))
            group_samples[group].add_file(FileInfo(sf.fullpath, FileContext.from_origin(sf.role+'_rep'+str(c))))
        
    
    #lets alleviate the load of generating pooled bams if they already exist.
    sys.stdout.write('Generating pooled samples.\n')
    sys.stdout.flush()
    group_samples_to_run = []
    group_samples_to_skip = []
    for s in group_samples.values():
        out_manifest_dest = os.path.join(s.dest, s.name+'_output_manifest.tsv')
        if os.path.isfile(out_manifest_dest):
            sys.stdout.write('-> It appears that %s has already been pooled. Skipping this sample...\n' % (s.name,))
            group_samples_to_skip.append(s)
        else:
            group_samples_to_run.append(s)
    sys.stdout.flush()
    
    run_read_manifest_pipeline(group_samples_to_skip)
    run_pipeline(make_replicate_pool_pipeline(args), group_samples_to_run, args)
    
    for gsample in (group_samples_to_run + group_samples_to_skip):
        gtemplate = None
        for s in samples:
            if s.name == gsample.get_attribute('template'):
                gtemplate = s
                break
        gsample_data = copy.deepcopy(gtemplate)
        gsample_data.set_attribute('pooledreplicate', True)
        gsample_data.name = gsample.name+'_pooled'
        for f in gsample.files:
            if f.cxt.module == 'Merge Treatment Bams':
                gsample_data.add_file(FileInfo(f, FileContext.from_origin('treatment')))
            elif f.cxt.module == 'Merge Control Bams':
                gsample_data.add_file(FileInfo(f, FileContext.from_origin('control')))
        samples.append(gsample_data)
    sys.stdout.write("Completed pooled replicate generation!\n\n")
    sys.stdout.write("--------------------------------------------\n\n")
    
    return samples
#end generate_replicate_pools()


def make_replicate_pool_pipeline(args):
    pool_pipeline = AnalysisPipeline('Generate Pooled Samples')
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferToShm
        x = TransferToShm.TransferToShm()
        x.set_parameter('shm_path', args.shm_path)
        pool_pipeline.append_module(x, critical=True)
        
    #merge bam files
    from ThackTech.Pipelines.PipelineModules import MergeBams
    x = MergeBams.MergeBams(name='Merge Treatment Bams', processors=args.threads)
    x.set_parameter('sort', True)
    x.set_parameter('index', True)
    x.set_parameter('postfix', 'treatment')
    x.set_resolver('alignments', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin() and f.context.role == 'treatment'))
    pool_pipeline.append_module(x, critical=True)
    
    
    x = MergeBams.MergeBams(name='Merge Control Bams', processors=args.threads)
    x.set_parameter('sort', True)
    x.set_parameter('index', True)
    x.set_parameter('postfix', 'control')
    x.set_resolver('alignments', lambda cxt: cxt.sample.find_files(lambda f: f.cxt.is_origin() and f.context.role == 'control'))
    pool_pipeline.append_module(x, critical=False)
    
    
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pool_pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
        
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    pool_pipeline.append_module(OutputManifest.OutputManifest())

    return pool_pipeline
#end make_replicate_pool_pipeline()


def generate_psuedoreplicate_samples(samples, args):
    sys.stdout.write('Generating sample pseudoreplicates.\n')
    sys.stdout.flush()
    
    psr_samples_to_run = []
    psr_samples_to_skip = []
    for s in samples:
        psr1_bam_dest = os.path.join(s.dest, 'pseudoreps', s.name+'.pr1.bam')
        psr2_bam_dest = os.path.join(s.dest, 'pseudoreps', s.name+'.pr2.bam')
        if os.path.isfile(psr1_bam_dest) and os.path.isfile(psr2_bam_dest):
            sys.stdout.write('-> It appears that %s already has pseudoreplicates generated. Skipping this sample...\n' % (s.name,))
            s.add_file('Pseudoreplicates', 'pr1', psr1_bam_dest)
            s.add_file('Pseudoreplicates', 'pr2', psr2_bam_dest)
            psr_samples_to_skip.append(s)
        else:
            psr_samples_to_run.append(s)
    sys.stdout.flush()
    
    run_read_manifest_pipeline(psr_samples_to_skip)
    run_pipeline(make_pseudoreplicate_pipeline(args), psr_samples_to_run, args)


    pseudo_reps = []
    for sample in (psr_samples_to_run + psr_samples_to_skip):
        #generate sample for psuedoreplicate 1
        psr1 = copy.deepcopy(sample)
        if psr1.has_attribute('pooledreplicate'):
            psr1.remove_attribute('pooledreplicate')
            psr1.set_attribute('pooledpseudoreplicate', True)
        else:
            psr1.set_attribute('pseudoreplicate', True)
        psr1.name = sample.name+'_psrep1'
        psr1.add_file('source', 'treatment', sample.get_file('Pseudoreplicates', 'pr1'))
        pseudo_reps.append(psr1)
        
        #generate sample for psuedoreplicate 2
        psr2 = copy.deepcopy(sample)
        if psr2.has_attribute('pooledreplicate'):
            psr2.remove_attribute('pooledreplicate')
            psr2.set_attribute('pooledpseudoreplicate', True)
        else:
            psr2.set_attribute('pseudoreplicate', True)
        psr2.name = sample.name+'_psrep2'
        psr2.add_file('source', 'treatment', sample.get_file('Pseudoreplicates', 'pr2'))
        pseudo_reps.append(psr2)
        
        #mark the sample as primary replicate if necessary (should be last to avoid dirtying the psr from copy)
        if not sample.has_attribute('pooledreplicate'):
            sample.set_attribute('primaryreplicate', True)
            
    samples = psr_samples_to_run + psr_samples_to_skip + pseudo_reps
    sample_count = len(samples)
    
    sys.stdout.write("Completed sample pseudoreplicates generation!\n\n")
    sys.stdout.write("--------------------------------------------\n\n")
    
#end generate_psuedoreplicate_samples()


def make_pseudoreplicate_pipeline(args):
    idr_pre_pipeline = AnalysisPipeline('Generate Pseudoreplicates')
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferToShm
        x = TransferToShm.TransferToShm()
        x.set_parameter('shm_path', args.shm_path)
        idr_pre_pipeline.append_module(x, critical=True)
        
    from ThackTech.Pipelines.PipelineModules import GeneratePseudoreplicates
    x = GeneratePseudoreplicates.GeneratePseudoreplicates()
    x.set_resolver('alignments', lambda s: s.get_file('source', 'treatment'))
    idr_pre_pipeline.append_module(x, critical=True)
    
    #sort pseudoreplicate files
    from ThackTech.Pipelines.PipelineModules import SortBam
    x = SortBam.SortBam()
    x.set_resolver('alignments', lambda s: s.get_file('Pseudoreplicates', 'pr1'))
    idr_pre_pipeline.append_module(x, critical=True)
    x = SortBam.SortBam()
    x.set_resolver('alignments', lambda s: s.get_file('Pseudoreplicates', 'pr2'))
    idr_pre_pipeline.append_module(x, critical=True)
    
    #index sorted merged bam file
    from ThackTech.Pipelines.PipelineModules import IndexBam
    x = IndexBam.IndexBam()
    x.set_resolver('alignments', lambda s: s.get_file('Pseudoreplicates', 'pr1'))
    idr_pre_pipeline.append_module(x, critical=True)
    x = IndexBam.IndexBam()
    x.set_resolver('alignments', lambda s: s.get_file('Pseudoreplicates', 'pr2'))
    idr_pre_pipeline.append_module(x, critical=True)
    
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        idr_pre_pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
        
    # from ThackTech.Pipelines.PipelineModules import OutputManifest
    # idr_pre_pipeline.append_module(OutputManifest.OutputManifest())

    return idr_pre_pipeline
#end make_pseudoreplicate_pipeline()


def run_idr_analysis(samples, args):
    idr_samples = []
    groups = set([s.get_attribute('group') for s in samples])
    for group in groups:
        samples_in_group = [s for s in samples if s.get_attribute('group') == group]
        gsample = PipelineSample(group, samples_in_group[0].genome, samples_in_group[0].dest)
        if samples_in_group[0].has_attribute('broad'):
            gsample.set_attribute('broad', True)
        for s in samples_in_group:
            if s.has_attribute('primaryreplicate'):
                gsample.add_file('primaryreplicate', s.name, resolve_idr_peak_file(s))
                
            elif s.has_attribute('pseudoreplicate'):
                gsample.add_file('pseudoreplicate', s.name, resolve_idr_peak_file(s))
                
            elif s.has_attribute('pooledpseudoreplicate'):
                gsample.add_file('pooledpseudoreplicate', s.name, resolve_idr_peak_file(s))
        idr_samples.append(gsample)

    for sample in idr_samples:
        sys.stderr.write(sample.name+'\n')
        sys.stderr.write(sample.files+'\n\n')
        
    sys.stdout.write('Performing IDR consistency analysis....\n')
    run_pipeline(make_IDR_analysis_pipeline(args), idr_samples, args)
            
    sys.stdout.write("Completed IDR consistency analysis!\n\n")
    sys.stdout.write("--------------------------------------------\n\n")
#end run_idr_analysis()


def make_IDR_analysis_pipeline(args):
    post_idr_pipeline = AnalysisPipeline('IDR Analysis')
    from ThackTech.Pipelines.PipelineModules import PerformIDRv2Analysis
    x = PerformIDRv2Analysis.PerformIDRv2Analysis()
    x.set_resolver('primary_replicates', lambda s: s.get_file_group('primaryreplicate').values())
    x.set_resolver('pseudo_replicates', lambda s: s.get_file_group('pseudoreplicate').values())
    x.set_resolver('pooled_pseudo_replicates', lambda s: s.get_file_group('pooledpseudoreplicate').values())
    x.set_parameter('primary_replicates_IDR_threshold', 0.01)
    x.set_parameter('pseudo_replicates_IDR_threshold', 0.01)
    x.set_parameter('pooled_pseudo_replicates_IDR_threshold', 0.0025)
    post_idr_pipeline.append_module(x, critical=True)
    
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    post_idr_pipeline.append_module(OutputManifest.OutputManifest())

    return post_idr_pipeline
#end make_IDR_analysis_pipeline()



def resolve_idr_peak_file(sample):
    if sample.has_file_group('MACS2'):
        if sample.has_attribute('broad') and sample.get_attribute('broad'):
            return sample.get_file('MACS2', 'broad_peaks')
        else:
            return sample.get_file('MACS2', 'narrow_peaks')
    else:
        return "" #right now only MACS2 peak calling is supported for IDR
#end resolve_idr_peak_file()

if __name__ == "__main__":
    main()

