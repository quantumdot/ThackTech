import os
import pandas as pd
import argparse
import sys
import copy
import multiprocessing
from ThackTech import Common
from ThackTech import filetools
from ThackTech import aligntools
from ThackTech.Pipelines import PipelineSample, AnalysisPipeline
from ThackTech.Pipelines.PipelineRunner import add_runner_args, get_configured_runner


#manifest should be structured as follows (with headers)
#Name    Group    Treatment    Control    Genome    Broad    Dest

gopts = {
    'shm_enabled': False,
    'conv_bdg_to_bw': True,
}

try:
    cpu_count = multiprocessing.cpu_count()
except:
    cpu_count = 1
if cpu_count is None:
    cpu_count = 1




class EpicPipelineSample(PipelineSample):

    def __init__(self, sample):
        PipelineSample.__init__(self, sample['Name'], sample['Genome'], sample['Dest'])
        self.add_file('source', 'treatment', sample['Treatment'])
        if 'Control' in sample:
            self.add_file('source', 'control', sample['Control'])
        if 'Broad' in sample and sample['Broad']:
            self.add_attribute('broad', True)
        if 'Group' in sample:
            self.add_attribute('group', sample['Group'])
    #end __init__()
#end PipelineSample




def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Given a sample manifest, runs Epic (SICER).',
        #epilog='The installed version of MACS is: '+show_macs_version(gopts['macs1_path'], None, False)+'; '+show_macs_version(gopts['macs2_path'], None, False)
    )
    parser.add_argument('manifest', help="Manifest file containing sample information. Manifest file should contain the following columns (tsv, case sensitive): [Name] [Genome] [Dest] [Treatment] [Control<optional>] [Broad<optional; bool-like>]")
    #parser.add_argument('--macs-version', choices=['macs1', 'macs2'], help="Specifies which version of MACS to run")
    parser.add_argument('-d', '--duplicates', default='auto', help="Specifies the MACS --keep-dup option. One of {'auto', 'all', <int>}.")
    parser.add_argument('-f', '--format', default=None, choices=["BED", "SAM", "BAM", "BAMPE", "BOWTIE", "ELAND", "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET"], help="Format of the alignment files for the entire run. If not specified, format will be auto-detected.")
    parser.add_argument('-bw', default=300, type=int, help="Bandwith (--bw) parameter for macs. Average sonnication fragment size expected from wet lab.")
    parser.add_argument('--sigout', default='bdg', choices=['wig', 'bdg'], help="Output type for signal. Either wig or bedgraph.")
    available_qc_choices = ['chance', 'fingerprint', 'frip']
    parser.add_argument('--qc', choices=available_qc_choices+['all'], default=[], action='append', help="Specify which QC routines to run. chance runs the IPStrength and spectrum modules from the chance package (Song et al. Genome Biology 2012. doi:10.1186/gb-2012-13-10-r98). fingerprint runs the DeepTools BamFingerprint module. FRiP computes the Fraction of Reads in Peaks. all will run all available QC modules.")
    parser.add_argument('--rpkmbw', action='store_true', help="Generate RPKM normailzed bigwig signal files.")
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
    
    samples = [EpicPipelineSample(s, args.format) for s in sample_manifest.to_dict(orient='records')]
    sample_count = len(samples)
    sys.stdout.write('-> Found %d item%s for processing.\n\n' % (sample_count, ('s' if sample_count > 1 else '')))
    
    
#end main()

def run_pipeline(pipeline, samples, args):
    runner = get_configured_runner(args, pipeline)
    runner.run(samples)
#end get_runner()


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
            x = MACS2Peakcall.MACS2Peakcall()
            x.set_parameter('duplicates', args.duplicates)
            x.set_parameter('bw', args.bw)
            pipeline.append_module(x, critical=True)
        else:
            from ThackTech.Pipelines.PipelineModules import MACS1Peakcall
            x = MACS1Peakcall.MACS1Peakcall()
            x.set_parameter('duplicates', args.duplicates)
            x.set_parameter('bw', args.bw)
            x.set_parameter('sigout', args.sigout)
            pipeline.append_module(x, critical=True)
            
        from ThackTech.Pipelines.PipelineModules import GenerateMACsModelFigure
        pipeline.append_module(GenerateMACsModelFigure.GenerateMACsModelFigure())
        
    from ThackTech.Pipelines.PipelineModules import GeneratePlainBed
    pipeline.append_module(GeneratePlainBed.GeneratePlainBed())
    
    if (args.sigout == 'bdg' or args.macs_version == 'macs2') and gopts['conv_bdg_to_bw']:
        from ThackTech.Pipelines.PipelineModules import ConvertBedgraphToBigWig
        pipeline.append_module(ConvertBedgraphToBigWig.ConvertBedgraphToBigWig())
        
    if 'chance' in args.qc:
        from ThackTech.Pipelines.PipelineModules import ChanceAnalysis
        pipeline.append_module(ChanceAnalysis.ChanceAnalysis())
        
    if 'fingerprint' in args.qc:
        from ThackTech.Pipelines.PipelineModules import BamFingerprint
        pipeline.append_module(BamFingerprint.BamFingerprint())
        
    if 'frip' in args.qc:
        from ThackTech.Pipelines.PipelineModules import FRiPAnalysis
        x = FRiPAnalysis.FRiPAnalysis()
        def resolve_frip_input(s):
            if s.has_file_group('MACS1'):
                return os.path.join(s.dest, s.name+'_peaks.bed')
            elif s.has_file_group('MACS2'):
                if s.has_attribute('broad') and s.get_attribute('broad'):
                    return os.path.join(s.dest, s.name+'_peaks.broadPeak')
                else:
                    return os.path.join(s.dest, s.name+'_peaks.narrowPeak')
            else:
                return None
        #end resolve_frip_input()
        x.set_resolver('peaks', resolve_frip_input)
        pipeline.append_module(x)
        
    if args.shm:
        from ThackTech.Pipelines.PipelineModules import TransferFromShm
        pipeline.append_module(TransferFromShm.TransferFromShm(), critical=True)
        
    if args.rpkmbw:
        from ThackTech.Pipelines.PipelineModules import BamToRpkmNormBigWig
        x = BamToRpkmNormBigWig.BamToRpkmNormBigWig()
        x.set_resolver('bam', lambda s: s.get_file('source', 'treatment'))
        pipeline.append_module(x)    
        
    from ThackTech.Pipelines.PipelineModules import OutputManifest
    pipeline.append_module(OutputManifest.OutputManifest())
    
    return pipeline
#end make_peak_calling_and_qc_pipeline




if __name__ == "__main__":
    main()

