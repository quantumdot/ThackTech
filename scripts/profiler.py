#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')

import os
import re
import sys
import time
import argparse
import metaseq
import metaseq.colormap_adjust as colormap_adjust
import pybedtools
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats
import multiprocessing
from ThackTech import filetools
from ThackTech.Plotting import sigcollector, stats as ttstats
from ThackTech import chromtools

#http://matplotlib.1069221.n5.nabble.com/Auto-wrapping-text-within-a-plot-Is-there-a-simpler-solution-td20.html



gopts = {
    'args': None,
    'co': None,
    'output_base': '',
    'fig_cols': 0,
    'fig_rows': 0,
    'x_axis': None,
    'plot_axes': {},
    'extra_artists': [],
}

plot_fonts = {
    'legend': None,
    'axis': None
}

class ProfileSample:
    def __init__(self, sample_id, sig_id, bed_id, signal_array, sig_label, bed_label):
        self.id = sample_id
        self.sig_id = sig_id
        self.bed_id = bed_id
        self.signal_array = signal_array
        self.sig_label = sig_label.replace("\\n", "\n")
        self.bed_label = bed_label.replace("\\n", "\n")
        self.show_yaxis = False
        self.hlines = []
    #end __init__()
#end class ProfileSample


def main():
    system_cpu_count = multiprocessing.cpu_count()

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    input_group = parser.add_argument_group('Input Data Options')
    input_group.add_argument('-b', '--bed', action='append', required='true', help='Bed files to profile against.')
    input_group.add_argument('-s', '--sig', action='append', required='true', help='Signal files to profile. BigWig, BAM, BED, BigBed files are supported, but all but bigwig will be treated by counting intervals.')
    input_group.add_argument('-i', '--inp', action='append', default=[], help='Input files to be used for normalization.')
    input_group.add_argument('--ignoreinput', action='store_true', help='Do not use input data.')
    
    
    labeling_group = parser.add_argument_group('Labeling Options')
    labeling_group.add_argument('-sl', '--slabel', action='append', default=[], help='Signal labels for each plot. Specify in the same order as the signal files. If not supplied, file basename will be used.')
    labeling_group.add_argument('-il', '--ilabel', action='append', default=[], help='Interval labels for each plot. Specify in the same order as the interval files. If not supplied, file basename will be used.')
    labeling_group.add_argument('--title', action='store', default='', help='Title for entire plot')
    labeling_group.add_argument('--bedscorelabel', action='store', default="Score", help='Label for plotting the score column of intervals.')
    labeling_group.add_argument('--xnumticks', action='store', type=int, default=3, help='Number of tick marks to use in the x-axis.')
    labeling_group.add_argument('--ynumticks', action='store', type=int, default=4, help='Number of tick marks to use in the y-axis.')
    labeling_group.add_argument('--xlabelrot', action='store', type=int, default=0, help='Angle, in degrees, for x axis label rotation.')
    labeling_group.add_argument('--ylabelrot', action='store', type=int, default=90, help='Angle, in degrees, for y axis label rotation.')
    
    profile_group = parser.add_argument_group('Profiling Options')
    profile_group.add_argument('--nochromfilter', action='store_true', help='Do not filter bed file for common chromosomes. By default profiling only occurs on chromosomes common to all files (useful for ignoring random or Un* chromosomes).')
    profile_group.add_argument('--chrignore', action='store', default=None, help='Regex for chromosome strings to filter out of intervals.')
    profile_group.add_argument('--up', action='store', default=1000, type=int, help='Span upstream (in bp) from selected intervals alignement to profile.')
    profile_group.add_argument('--down', action='store', default=1000, type=int, help='Span downstream (in bp) from selected intervals alignement to profile.')
    profile_group.add_argument('--res', action='store', default=10, type=int, help='Profiling resolution in bp.')
    profile_group.add_argument('--dir', action='store_true', help='If set, the direction (+/-) [strand orientation] is considered in profiling. Strand information should be present in bed file used, and if not present is assumed to be +.')
    profile_group.add_argument('--align', action='store', choices=['left', 'center', 'right', 'scale'], default='center', help='Method used to align intervals. Left aligns intervals along the 5\' end, while right aligns intervals aling the 3\' end (assuming intervals provide strand information and --dir is specified, otherwise + strand is assumed). Center aligns intervals along the center center point of intervals. Scale will scale all intervals to the same apparent width specified by --scaleregionsize.')
    profile_group.add_argument('--scaleregionsize', action='store', type=int, default='3000', help='Size of the scaled region (i.e. gene-body) in bp. This option is only useful for "--align scale" option.')
    profile_group.add_argument('--nan', action='store', choices=['zero', 'ignore'], default='ignore', help='How to handle missing or NaN values.')
    profile_group.add_argument('--collectionmethod', action='store', choices=['get_as_array', 'ucsc_summarize', 'summarize'], default='get_as_array', help='Method for collecting signal data.')
    
    scale_group = parser.add_argument_group('Scaling Options')
    scale_group.add_argument('--scalegroups', action='store', default=None, help='Groups of plots to share color/y-axis scales. If not specified, all plots will be constructed with the same scale. Parameter should be specified as comma-separated lists of 0-based offsets of samples, and groups separated with a semicolon. Ex: 0;1,2;3,4 results in sample 0 plotted with independent scale, 1 and 2 sharing scale, and 3 and 4 sharing scale. If specified, but parameter omits samples, then the omitted samples will each be scaled independently.')
    scale_group.add_argument('--scalebedgroups', action='store', default=None, help='Groups of plots to share color/y-axis scales. If not specified, all plots will be constructed with the same scale. Parameter should be specified as comma-separated lists of 0-based offsets of samples, and groups separated with a semicolon. Ex: 0;1,2;3,4 results in sample 0 plotted with independent scale, 1 and 2 sharing scale, and 3 and 4 sharing scale. If specified, but parameter omits samples, then the omitted samples will each be scaled independently.')
    scale_group.add_argument('--saturatemin', action='store', type=float, default=0.01, help='In the heatmap plot, saturate the <--saturatemin> percent bottom values.')
    scale_group.add_argument('--saturatemax', action='store', type=float, default=0.01, help='In the heatmap plot, saturate the <--saturatemax> percent top values.')
    scale_group.add_argument('--coefficients', action='store', nargs='*', type=float, help='Coefficients to multiply signals by, one value for each signal submitted.')
    scale_group.add_argument('--normalizationmethod', action='store', default='log2', choices=['ratio', 'log2', 'reciprocal_ratio', 'subtract', 'add', 'mean'], help='Method to use for normalizing signal by input.')
    
    clustsort_group = parser.add_argument_group('Clustering/Sorting Options')
    clustsort_group.add_argument('--sort', action='store', default=None, type=int, help='sample index (0-based) to use for sorting. If not specified than the order of the bed file is used. Mutually exclusive with --kmeans.')
    clustsort_group.add_argument('--sortmethod', action='store', choices=['mean', 'median', 'max', 'min', 'sum'], default='mean', help='Method used for sorting.')
    clustsort_group.add_argument('--sortrange', action='store', default=None, help='Range of the profiles (in relative bp) to be used in the sorting operation. Specify in the format "start:stop". Default is to use the entire range.')
    clustsort_group.add_argument('--kmeans', action='store_true', help='If set, Perform K-means clustering on the data. Mutually exclusive with --sort.')
    clustsort_group.add_argument('--k', action='store', type=int, help='Number of clusters, k, to fit data to when performing K-means clustering.')
    clustsort_group.add_argument('--autok', action='store_true', help='Optimize the number of clusters, k, in the dataset. Mutually exclusixe with --k.')
    clustsort_group.add_argument('--ksamples', action='store', default='all', help='Comma-separated list of 0-based offsets of samples to use for K-means clustering. Use \'all\' to cluster on all samples.')
    clustsort_group.add_argument('--summarymethod', action='store', choices=['mean', 'median', 'max', 'min', 'sum'], default='mean', help='Method used for producing summary (avg) plot.')
    clustsort_group.add_argument('--hline', action='store_true', help='Draw horizontal lines to delineate clusters.')
    clustsort_group.add_argument('--hlineweight', action='store', type=float, default=0.5, help='Horizontal line weight.')
    clustsort_group.add_argument('--hlinestyle', action='store', default="-", help='Line style for the horizontal lines. See matplotlib axhline documentation for valid values.')
    
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--name', action='store', default='', help='Base name for the plot output.')
    output_group.add_argument('--format', action='store', default='pdf', choices=['pdf', 'png', 'svg'], help='Format to output the final figure.')
    output_group.add_argument('--plot', action='append', choices=['avg', 'violin','heat','kavg','bedscores','truebedscores'], required='true', help='Types of plots to produce. Supply multiple options to produce hybrid plots. avg will produce average profiles, heat will produce heatmaps, kavg will produce average profiles for each class determined by k-means clustering (only available when used with --plot heat AND --kmeans). bedscores will plot the score column for intervals simply, while truebedscores will plot the score column for intervals more accuratly (showing true genomic context of element; mutually exclusive with --plot bedscores).')
    output_group.add_argument('--dpi', action='store', type=int, default=600, help='DPI resolution of the saved figure.')
    output_group.add_argument('--width', action='store', type=float, default=8.0, help='Width (in inches) of the figure.')
    output_group.add_argument('--height', action='store', type=float, default=6.0, help='Height (in inches) of the figure.')
    output_group.add_argument('--avgplotrows', action='store', type=int, default=1, help='Number of rows to use for average profile plots.')
    output_group.add_argument('--heatplotrows', action='store', type=int, default=3, help='Number of rows to use for heatmap plots.')
    output_group.add_argument('--rotate', action='store_true', help='By default plots will be arranged with signals going across and intervals going down. --rotate will change the orientation so that intervals go across and signals going down.')
    output_group.add_argument('--dumpsummary', action='store_true', help='Write summary data to a text file. See --summaryrange and --summarymethod for details on affecting the output generated by this switch.')
    output_group.add_argument('--summaryrange', action='store', default=None, help='Range(s) of the profiles (in relative bp) to be used in the dump summary operation. Specify in the format "start:stop". Default is to use the entire range. Write multiple ranges by separating ranges with a semicolon')
    #output_group.add_argument('--summarymethod', action='store', choices=['mean', 'median', 'max', 'min', 'sum'], default='mean', help='Method used for producing summary data.')
    
    
    plot_group = parser.add_argument_group('General Plotting Options')
    plot_group.add_argument('--vline', action='store_true', help='Draw a vertical line at the "zero" point.')
    plot_group.add_argument('--vlineweight', action='store', type=float, default=0.5, help='Vertical line weight.')
    plot_group.add_argument('--vlinestyle', action='store', default=":", help='Line style for the vertical line. See matplotlib axvline documentation for valid values.')
    plot_group.add_argument('--fontsize', action='store', type=int, default=8, help='Font size used in plots.')
    plot_group.add_argument('--legendfontsize', action='store', type=int, default=6, help='Font size used in plot legend.')
    plot_group.add_argument('--showci', action='store_true', help='Plot confidence interval for average plot.')
    plot_group.add_argument('--ciwidth', action='store', default='sem', help='Confidence interval width for average plot. sem for Standard Error of the mean, std for Standard deviation, or float for percent CI (i.e. 0.95).')
    plot_group.add_argument('--genome', action='store', help='UCSC reference genome that intervals and signal are computed with (i.e. mm9, hg19, etc.).')
    plot_group.add_argument('--colors', action='store', nargs='+', default=['r','g','b','c','m','y','k','firebrick','darkolivegreen','navy','palevioletred'], help="Colors to cycle through for individual samples on aggregate plots.")
    plot_group.add_argument('--heatcolor', action='store', default='#b11902', help="Color of the max positive heatmap value, used to derive heatmap colormap that is centered on zero, and accounts for asymmetrical vmin and vmax by matching saturation/value of high and low colors.")
    plot_group.add_argument('--linewidth', action='store', type=float, default=0.5, help="Width of line used for aggregate plots (i.e. average profile).")
    
    resources_group = parser.add_argument_group('Resource Options')
    resources_group.add_argument('--cache', action='store_true', help='If set, profiles are dumped to disk.')
    resources_group.add_argument('--cachedir', action='store', default=".profiler_cache", help='Location to place cache files.')
    resources_group.add_argument('-p', '--processors', action='store', dest='cpus', metavar='cpus', default=(system_cpu_count/2), type=int, help='Number of processors to use.')
    gopts['args'] = args = parser.parse_args()
    document_args()

    if len(args.slabel) == 0:
        sys.stderr.write("WARNING: Signal labels were not supplied. Using signal file basenames instead.\n")
    if len(args.ilabel) == 0:
        sys.stderr.write("WARNING: Interval labels were not supplied. Using interval file basenames instead.\n")
    if len(args.slabel) > 0 and len(args.slabel) < len(args.sig):
        sys.stderr.write("ERROR: Fewer signal labels supplied than there are signals! Signal and label counts must be equal if supplying labels!\n")
        sys.exit(1)
    if len(args.ilabel) > 0 and len(args.ilabel) < len(args.bed):
        sys.stderr.write("ERROR: Fewer interval labels supplied than there are intervals! Interval and label counts must be equal if supplying labels!\n")
        sys.exit(1)
    if len(args.inp) > 0 and len(args.inp) < len(args.sig):
        sys.stderr.write("ERROR: Fewer inputs supplied than there are samples! sample and input counts must be equal if supplying inputs!\n")
        sys.exit(1)
    if len(args.plot) <= 0:
        sys.stderr.write("ERROR: No plot types were selected!\n")
        sys.exit(1)
    if args.cpus <= 0:
        sys.stderr.write("WARNING: Cannot use fewer than one processor for computation. Setting --processors = 1!\n")
        args.cpus = 1
    if args.cpus > system_cpu_count:
        sys.stderr.write("WARNING: Cannot use more than the number of processors on the system ("+system_cpu_count+") for computation. Setting --processors = "+system_cpu_count+"!\n")
        args.cpus = system_cpu_count
    if len(args.bed) > 1 and args.kmeans:
        sys.stderr.write("ERROR: Cannot perform K-means analysis when using multiple beds!\n")
        sys.exit(1)
    if 'kavg' in args.plot and not args.kmeans and not 'heat' in args.plot:
        sys.stderr.write("WARNING: to plot averages for clusters, --kmeans and --plot heat must be used! Ignoring this directive!\n")
        args.plot.remove('kavg')
    
    
    
    collection_opts = sigcollector.CollectorOptions()
    collection_opts.align = args.align
    collection_opts.upstream = abs(args.up)
    collection_opts.downstream = abs(args.down)
    collection_opts.scaleregionsize = abs(args.scaleregionsize)
    collection_opts.resolution = abs(args.res)
    collection_opts.direction = bool(args.dir)
    collection_opts.validate()
    gopts['co'] = collection_opts
    gopts['x_axis'] = collection_opts.xaxis
    
    gopts['output_base'] = "%s.%s" % (args.name, str(collection_opts))
    gopts['savename_notes'] = []
    
    
    if args.ignoreinput:
        gopts['args'].inp = args.inp = []
    
    
    
    gopts['chromsets'] = chromtools.get_common_chroms(args.bed, args.sig)
    if not args.nochromfilter:
        sys.stderr.write("Filtering common chromomsomes....\n")
        gopts['chromsets'].use = gopts['chromsets'].common
        sys.stderr.write("=>  Keeping: "+', '.join(gopts['chromsets'].common)+"\n")
        sys.stderr.write("=> Ignoring: "+', '.join(gopts['chromsets'].uncommon)+"\n")
        sys.stderr.write("\n")
    else:
        gopts['chromsets'].use = gopts['chromsets'].all
    if len(gopts['chromsets'].use) < 1:
        sys.stderr.write("ERROR: No valid chromosomes left after filtering for only common chromosomes!\n")
        sys.stderr.write(" -> Try using the --nochromfilter or --chrignore options.\n")
        sys.exit(1)
        
    if args.chrignore is not None:
        pattern = re.compile(args.chrignore)
        filtered_chrs = set()
        for chrom in gopts['chromsets'].use:
            if pattern.search(chrom) is not None:
                filtered_chrs.add(chrom)
        sys.stderr.write('Filtering chromomsomes matching filter "%s"....\n' % (args.chrignore,))
        sys.stderr.write("=>  Ignoring: "+', '.join(filtered_chrs)+"\n")
        gopts['chromsets'].use -= filtered_chrs
        sys.stderr.write("\n")
        
        
    
    samples = []
    collection_args = {
        'norm_method': args.normalizationmethod,
        'norm_pseudocount': 1.0,
        'cache_dir': (args.cachedir if args.cache else None), 
        'cache_base': args.name, 
        'collectionmethod': args.collectionmethod, 
        'cpus': args.cpus 
    }
    for s in xrange(len(args.sig)):
        for b in xrange(len(args.bed)):
            sys.stderr.write("Processing %s vs %s....\n" % (args.bed[b], args.sig[s]))
            sys.stderr.write("-> Preparing intervals.....\n")
            
            bedtool = sigcollector.IntervalProvider(args.bed[b], collection_opts, args.genome, gopts['chromsets'].use)
            #bedtool1 = expand_bed(args.up, args.down, args.align, , gopts['chromsets'].use)
            #bedtool2 = expand_bed(args.up, args.down, args.align, args.bed[b], gopts['chromsets'].use)#seems to choke using the same instance...... not sure why!!!
            s_label = args.slabel[s] if len(args.slabel)-1 >= s else os.path.splitext(os.path.basename(args.sig[s]))[0]
            b_label = args.ilabel[b] if len(args.ilabel)-1 >= b else os.path.splitext(os.path.basename(args.bed[b]))[0]
            input_sig = args.inp[s] if len(args.inp)-1 >= s else None
            if input_sig is not None and input_sig.lower() == 'none':
                input_sig = None
            
            signal = sigcollector.get_signal(bedtool, s_label+b_label, args.sig[s], input_sig, **collection_args)
            
            ps = ProfileSample(len(samples), s, b, signal, s_label, b_label)
            ps.bedtool = bedtool
            #sys.stderr.write("NaN count: %d\n" % (np.isnan(ps.signal_array).sum(),))
            ps.signal_array = ttstats.correct_invalid_data(ps.signal_array, args.nan)
            #sys.stderr.write("NaN count: %d\n" % (np.isnan(ps.signal_array).sum(),))
            samples.append(ps)
            sys.stderr.write("\n")
    
    
    if 'bedscores' in args.plot or 'truebedscores' in args.plot:
        for b in xrange(len(args.bed)):
            bedtool = sigcollector.IntervalProvider(args.bed[b], collection_opts, args.genome, gopts['chromsets'].use)
            if 'truebedscores' in args.plot:
                signal = sigcollector.get_bed_score_signal_complex(bedtool, **collection_args)
            else:
                signal = sigcollector.get_bed_score_signal(bedtool)
            #print signal
            b_label = args.ilabel[b] if len(args.ilabel)-1 >= b else os.path.splitext(os.path.basename(args.bed[b]))[0]
            ps = ProfileSample(len(samples), len(args.sig), b, signal, args.bedscorelabel, b_label)
            ps.bedtool = bedtool
            samples.append(ps)
        gopts['savename_notes'].append("bed-score")
    
    if args.coefficients is not None:
        for s in samples:
            c = 1.0
            if len(args.coefficients) > s.sig_id:
                c = args.coefficients[s.sig_id]
            sys.stderr.write("Scaling sample {} by factor {}\n".format(s.sig_id, c))
            s.signal_array *= c
    
    
    gopts['group_count'] = count_groups(args.scalegroups, samples)
    if args.rotate:
        gopts['fig_cols'] = len(args.bed) + 1#gopts['group_count']
        gopts['fig_rows'] = 0
        if 'heat' in args.plot:
            gopts['fig_rows'] += (args.heatplotrows * len(args.sig))
        if 'avg' in args.plot:
            gopts['fig_rows'] += args.avgplotrows
        if 'violin' in args.plot:
            gopts['fig_rows'] += args.avgplotrows
    else:
        gopts['fig_cols'] = len(args.sig) + gopts['group_count']
        gopts['fig_rows'] = 0
        if 'heat' in args.plot:
            gopts['fig_rows'] += (args.heatplotrows * len(args.bed))
        if 'avg' in args.plot:
            gopts['fig_rows'] += args.avgplotrows
        if 'violin' in args.plot:
            gopts['fig_rows'] += args.avgplotrows
    
    
    if 'bedscores' in args.plot or 'truebedscores' in args.plot:
        if args.rotate:
            gopts['fig_rows'] += args.heatplotrows
        else:
            gopts['fig_cols'] += 1 #allocation for cbar was made by count_groups
        
        
    
    plt.rcParams['font.size'] = args.fontsize
    plt.rcParams['legend.fontsize'] = args.legendfontsize
    plt.rcParams['pdf.fonttype'] = 42 # use TrueType fonts
    plt.rcParams['ps.fonttype'] = 42
    fig = plt.figure(figsize=(args.width, args.height), dpi=args.dpi)
    
    
    #compute sort orders/grouping/clustering
    if args.kmeans:
        if args.ksamples == 'all':
            kindicies = range(len(samples))
        else:
            kindicies = [int(i) for i in args.ksamples.split(',')]
        if args.autok:
            k = range(1,10)
        else:
            k = args.k
        compute_Kmeans(samples, kindicies, k, args.bed[0], gopts['chromsets'].use)
        gopts['savename_notes'].append("k%d" %(gopts['k_info']['k'],))
    else:
        if args.sortrange is not None:
            start, stop = args.sortrange.split(':')
            start = int(int(start) / args.res)
            stop = int(int(stop) / args.res)
        else:
            start = 0
            stop = gopts['co'].total_bins
        compute_sorting(samples, args.sort, args.sortmethod, (start, stop))
        if args.sort is None:
            gopts['savename_notes'].append("sort-bed")
        else:
            gopts['savename_notes'].append("sort%d%s" % (args.sort, args.sortmethod))
        
    #compute saturation points and average profile limits
    compute_group_scales(args.scalegroups, samples, args.saturatemin, args.saturatemax)
    
    #generate the subplots....
    sys.stderr.write("Generating Figure....\n")
    for s in samples:
        add_signal_to_figure(s)
    sys.stderr.write("\n")
    
    #if we have multiple bed and multiple signals, add legend outside last avg plot
    if 'avg' in gopts['args'].plot and (len(args.sig) > 1 or len(args.bed) > 1):
        last_avg_ax = get_plot_axes('leg', 0, 0, 0)
        leg = last_avg_ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        leg.get_frame().set_linewidth(0.1)
        gopts['extra_artists'].append(leg)
    
    #add color scale bar if necessary
    if 'heat' in gopts['args'].plot:
        make_colormap_strip_for_groups(fig, args.scalegroups, samples)
    
    #add sensible x-axis label
    gopts['extra_artists'].append(fig.text(0.5, 0.04, collection_opts.xaxis_label, ha='center', va='center'))
    gopts['extra_artists'].append(fig.suptitle(args.title))
        
    #finally save the figure!
    save_figure(fig, "_".join(gopts['savename_notes']))
    plt.close(fig)
    
    if args.dumpsummary:
        ranges = []
        if args.summaryrange is not None:
            sr = args.summaryrange.split(';')
            for r in sr:
                start, stop = r.split(':')
                ranges.append((int(int(start) / args.res), int(int(stop) / args.res)))
        else:
            ranges.append((0, gopts['co'].total_bins))
        write_summary_profiles(samples, ranges, args.summarymethod)
    sys.stderr.write('Done!')
#end main()

def document_args():
    
    sys.stderr.write("\n")
    sys.stderr.write("# ARGUMENTS:\n")
    sys.stderr.write("# Run Start: %s\n" % (time.strftime("%a, %d %b %Y %H:%M:%S"),))
    sys.stderr.write("# Name: %s\n" % (gopts['args'].name,))
    
    sys.stderr.write("# Interval File(s):\n")
    for i in range(len(gopts['args'].bed)):
        sys.stderr.write("#\t%d: %s\n" % (i, os.path.abspath(gopts['args'].bed[i]),))
        
    sys.stderr.write("# Signal File(s):\n")
    for i in range(len(gopts['args'].sig)):
        sys.stderr.write("#\t%d: %s\n" % (i, os.path.abspath(gopts['args'].sig[i]),))
    
    sys.stderr.write("# Input File(s):\n")
    if len(gopts['args'].inp) > 0:
        for i in range(len(gopts['args'].inp)):
            if gopts['args'].inp[i].lower() == 'none':
                sys.stderr.write("#\t%d: None\n" % (i,))
            else:
                sys.stderr.write("#\t%d: %s\n" % (i, os.path.abspath(gopts['args'].inp[i]),))
    else:
        sys.stderr.write("#\tNone\n")
        
    sys.stderr.write("# --------------------------------------------------\n")
    sys.stderr.write("# Alignment: %s\n" % (gopts['args'].align,))
    if gopts['args'].align == 'scale':
        sys.stderr.write("# Span: %dbp Up; %dbp Meta; %dbp Down\n" % (gopts['args'].up, gopts['args'].scaleregionsize, gopts['args'].down))
    else:
        sys.stderr.write("# Span: %dbp Up; %dbp Down\n" % (gopts['args'].up, gopts['args'].down))
    sys.stderr.write("# Resolution: %d bp\n" % (gopts['args'].res,))
    sys.stderr.write("# Direction: %s\n" % ('ON' if gopts['args'].dir else 'OFF',))
    sys.stderr.write("# Filter common chromosomes: %s\n" % ('OFF' if gopts['args'].nochromfilter else 'ON',))
    sys.stderr.write("# Scale Groups: %s\n" % (gopts['args'].scalegroups,))
    if gopts['args'].kmeans:
        sys.stderr.write("# Cluster by K-means:\n")
        if gopts['args'].autok:
            sys.stderr.write("#    -> K = auto\n")
        else:
            sys.stderr.write("#    -> K = %d\n" % (gopts['args'].k,))
        sys.stderr.write("#    -> Cluster on samples: %s\n" % (gopts['args'].ksamples,))
    else:
        if gopts['args'].sort is None:
            sys.stderr.write("# Sort By: BED order\n")
        else:
            sys.stderr.write("# Sort By: %s of signal #%d\n" % (gopts['args'].sortmethod, gopts['args'].sort,))
    sys.stderr.write("# Saturate Min: %f\n" % (gopts['args'].saturatemin,))
    sys.stderr.write("# Saturate Max: %f\n" % (gopts['args'].saturatemax,))
    sys.stderr.write("# Plot Types: %s\n" % (", ".join(gopts['args'].plot),))
    sys.stderr.write("# Cache Data: %s\n" % ('ON' if gopts['args'].cache else 'OFF',))
    sys.stderr.write("# Num Threads: %d\n" % (gopts['args'].cpus,))
    sys.stderr.write("\n")
#end document_args()



def compute_Kmeans(samples, kindicies, k, bed, white_chroms=None):
    sys.stderr.write("Computing K-means clustering....\n")
    sys.stderr.write("-> Considering samples %s (0-based)\n" % (str(kindicies),))
    sys.stderr.write("-> Using K = %s\n" % (str(k),))
    
    pooled = np.hstack(tuple([s.signal_array for s in samples if s.id in kindicies]))
    ind, breaks = metaseq.plotutils.new_clustered_sortind(pooled, k=k, row_key=np.mean, cluster_key=np.median)
    
    if not isinstance(k, int):
        sys.stderr.write("    => Using auto K - found %d clusters\n" % (len(breaks),))
    
    hlines = []
    for b in breaks:
        hlines.append(b)
        
    for s in samples:
        s.hlines = hlines
        s.sort_order = np.argsort(ind)
    
    classes = make_interval_classes(ind, breaks, bed, white_chroms)
    for i in range(len(breaks)):
        sys.stderr.write("    -> Cluster %d: %d\n" % ((i+1), len(classes[(classes['class_id'] == (i+1))])))
    
    sys.stderr.write('\n')
    gopts['k_info'] = {
        'k': len(breaks),
        'ind': ind,
        'breaks': breaks,
        'classes': classes
    }
    return classes
#end compute_Kmeans

def make_interval_classes(sort_indicies, breaks, bed, white_chroms=None):
    intervals = pd.read_csv(bed, sep='\t', comment='#', skip_blank_lines=True, header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    if white_chroms is not None:
        intervals.drop(intervals[~intervals['chrom'].isin(white_chroms)].index, inplace=True)
    intervals.reset_index(inplace=True)
    intervals['sort_order'] = sort_indicies
    intervals.sort_values(by='sort_order', ascending=True, inplace=True)
    total_intervals = len(intervals)
    k = len(breaks)
    #print breaks
    classes = []
    for i in reversed(range(k-1)):
        size = breaks[i+1] - breaks[i]
        classes.extend([k-(i+1)] * size)
    classes.extend([k] * breaks[0])
    intervals['class_id'] = classes
    
    intervals.to_csv(gopts['output_base']+'.kmeans_classes.tsv', sep='\t', index=False)
    intervals.sort_index(inplace=True)
    return intervals
#end make_interval_classes()

def write_summary_profiles(samples, ranges, method):
    iv_beds = {}
    for s in samples:
        if s.bed_label not in iv_beds:
            iv_beds[s.bed_label] = pd.DataFrame(s.bedtool.generate_origional_as_tuples(), columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    for s in samples:
        b = iv_beds[s.bed_label]
        for r in ranges:
            b["{}_[{}..{}]".format(s.sig_label, r[0], r[1])] = ttstats.summarize_data(s.signal_array, method=method, srange=r, axis=1)
    
    for b in iv_beds:
        iv_beds[b].to_csv("{}.summary_{}.tsv".format(gopts['output_base'], b), sep='\t', index=False)

#end make_interval_classes()

def compute_sorting(samples, sort_index, method, sort_range):
    #print range
    if sort_index is None:
        for s in samples:
            s.sort_order = None
    else:
        sort_orders = {} #sort orders indexed by interval id
        for s in samples:
            if s.sig_id == sort_index:
                sort_orders[s.bed_id] = ttstats.summarize_data(s.signal_array, method=method, srange=sort_range, axis=1)
        for s in samples:
            s.sort_order = sort_orders[s.bed_id]
#end compute_sorting()




def get_groups(groups, samples):
    if groups is None:
        return [list(set(str(s.sig_id) for s in samples))]
    else:
        non_covered_samples = set([str(s.sig_id) for s in samples]) - set(groups.replace(";",",").split(","))
        if len(non_covered_samples) > 0:
            groups += ";" + ";".join(non_covered_samples)
        final_groups = []
        for group in groups.split(";"):
            final_groups.append([str(s.sig_id) for s in samples if str(s.sig_id) in group.split(",")])
        return [sl for sl in final_groups if len(sl) > 0]
#end get_groups()

def count_groups(groups, samples):
    return len(get_groups(groups, samples))
#end count_groups()

def compute_group_scales(groups, samples, min_saturate, max_saturate):
    groups_list = get_groups(groups, samples)
    #print groups_list
    for i in range(len(groups_list)):
        compute_scales_for_group(i, [s for s in samples if str(s.sig_id) in groups_list[i]], min_saturate, max_saturate)
#end compute_group_scales()

def compute_scales_for_group(group_id, samples, min_saturate, max_saturate):
    pooled = np.vstack(tuple([s.signal_array for s in samples]))
    heat_min = np.percentile(pooled[np.isfinite(pooled)], (min_saturate*100))
    heat_max = np.percentile(pooled[np.isfinite(pooled)], ((1-max_saturate)*100))
    avg_min = float('inf')
    avg_max = float('-inf')
    raw_min = float('inf')
    raw_max = float('-inf')
    for s in samples:
        mean_array = ttstats.summarize_data(s.signal_array, method=gopts['args'].summarymethod, axis=0)
        summary_array = ttstats.summarize_data(s.signal_array, method=gopts['args'].summarymethod, axis=1)
        t_mean_min = mean_array.min()
        t_mean_max = mean_array.max()
        t_raw_min = summary_array.min()
        t_raw_max  = summary_array.max()
        if t_mean_min < avg_min:
            avg_min = t_mean_min
        if t_mean_max > avg_max:
            avg_max = t_mean_max
            
        if t_raw_min < raw_min:
            raw_min = t_raw_min
        if t_raw_max > raw_max:
            raw_max = t_raw_max
    
    avg_max = avg_max + (abs(avg_max) * 0.1)
    if avg_min >= 0:
        avg_min = 0
    else:
        avg_min = avg_min - (abs(avg_min) * 0.1)
        
    raw_max = raw_max + (abs(raw_max) * 0.1)
    if raw_min >= 0:
        raw_min = 0
    else:
        raw_min = raw_min - (abs(raw_min) * 0.1)
    
    for s in samples:
        s.group = group_id
        if (gopts['args'].rotate and s.bed_id == samples[0].bed_id) or (not gopts['args'].rotate and s.sig_id == samples[0].sig_id):
            s.show_yaxis = True
        s.heat_min = heat_min
        s.heat_max = heat_max
        s.avg_min = avg_min
        s.avg_max = avg_max
        s.raw_min = raw_min
        s.raw_max = raw_max
        
    sys.stderr.write("-> Computed group %d heatmap min/max (%f, %f)\n" % (group_id, heat_min, heat_max))
    sys.stderr.write("-> Computed group %d average profile min/max (%f, %f)\n" % (group_id, avg_min, avg_max))
#end compute_scales_for_group()

def get_plot_axes(plot_type, group, bed_id, sig_id):
    #print "getting axis for %s (%d, %d, %d)\n" % (plot_type, group, bed_id, sig_id)
    if gopts['args'].rotate:
        col = bed_id #+ group
    else:
        col = sig_id + group
    
    if plot_type == 'heat':
        rowspan = gopts['args'].heatplotrows
        if gopts['args'].rotate:
            row = sig_id * gopts['args'].heatplotrows
        else:
            row = bed_id * gopts['args'].heatplotrows
    
    elif plot_type == 'avg':
        rowspan = gopts['args'].avgplotrows
        if 'violin' in gopts['args'].plot:
            row = gopts['fig_rows'] - (2 * gopts['args'].avgplotrows)
        else:
            row = gopts['fig_rows'] - gopts['args'].avgplotrows
        
    elif plot_type == 'violin':
        rowspan = gopts['args'].avgplotrows
        row = gopts['fig_rows'] - gopts['args'].avgplotrows
    
    elif plot_type == 'cbar':
        #special case: bed_id = min(sig_id) and sig_id = max(sig_id)
        if gopts['args'].rotate:
            rowspan = (sig_id - bed_id + 1) * gopts['args'].heatplotrows
            row = bed_id * gopts['args'].heatplotrows
            col = gopts['fig_cols']-1
        else:
            if 'violin' in gopts['args'].plot and 'avg' in gopts['args'].plot:
                rowspan = gopts['fig_rows'] - (2 * gopts['args'].avgplotrows)
            elif 'violin' in gopts['args'].plot or 'avg' in gopts['args'].plot:
                rowspan = gopts['fig_rows'] - gopts['args'].avgplotrows
            else:
                rowspan = gopts['fig_rows']
            row = 0
            col = sig_id + group + 1
    
    elif plot_type == 'leg':
        #we are tryign to return the bottom-right-most average plot
        #return get_plot_axes('avg', )
        rowspan = gopts['args'].avgplotrows
        
        if 'violin' in gopts['args'].plot and 'avg' in gopts['args'].plot:
            row = gopts['fig_rows'] - (2 * gopts['args'].avgplotrows)
            
        elif 'violin' in gopts['args'].plot or 'avg' in gopts['args'].plot:
            row = gopts['fig_rows'] - gopts['args'].avgplotrows
            
        else:
            row = gopts['fig_rows'] #should probably not be able to get here....
            
        if 'heat' in gopts['args'].plot:
            col = gopts['fig_cols'] - 2
        else:
            col = gopts['fig_cols'] - 1
            
            
    if (row,col) not in gopts['plot_axes']:
        ax = plt.subplot2grid((gopts['fig_rows'],gopts['fig_cols']), (row, col), rowspan)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=gopts['args'].xnumticks-1))
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=gopts['args'].ynumticks-1))
        [i.set_linewidth(0.1) for i in ax.spines.itervalues()]
        gopts['plot_axes'][(row,col)] = ax
    
    sys.stderr.write("getting axis for {} at ({}, {}) in grid ({}, {})\n".format(plot_type, row, col, gopts['fig_rows'],gopts['fig_cols']))
    return gopts['plot_axes'][(row,col)]
#end get_plot_axes()

def add_signal_to_figure(sample):
    #avg_profile_row = gopts['fig_rows']-gopts['args'].avgplotrows #3 if 'heat' in gopts['args'].plot else 0
    sys.stderr.write("-> Plotting data (%d, %d, %d) for %s x %s....\n" % (sample.group, sample.bed_id, sample.sig_id, sample.bed_label, sample.sig_label))
    
    if 'heat' in gopts['args'].plot:
        sys.stderr.write("    -> Generating heatmap...\n")
        #sys.stderr.write("    -> (%d, %d)\n" % (sample.bed_id*3, sample.sig_id))
        ax =  get_plot_axes('heat', sample.group, sample.bed_id, sample.sig_id)
        make_sig_heatmap(ax, sample)
        
    if 'violin' in gopts['args'].plot:
        sys.stderr.write("    -> Generating violin plot...\n")
        ax =  get_plot_axes('violin', sample.group, sample.bed_id, sample.sig_id)
        
        if gopts['args'].rotate:
            color = gopts['args'].colors[sample.sig_id % len(gopts['args'].colors)]
        else:
            color = gopts['args'].colors[sample.bed_id % len(gopts['args'].colors)]
        
        make_violin_plot(ax, sample, color)
        
    if 'avg' in gopts['args'].plot:
        ax =  get_plot_axes('avg', sample.group, sample.bed_id, sample.sig_id)
        if 'kavg' in gopts['args'].plot:
            sys.stderr.write("    -> Generating average profile for clusters...\n")
            #sys.stderr.write("    -> (%d, %d)\n" % (avg_profile_row,sample.sig_id))
            clust_colors = ['r','g','b','c','m','y']
            for i in range(gopts['k_info']['k']):
                sys.stderr.write("        -> Cluster #%d...\n" % ((i+1),))
                add_masked_group_to_avg_plot(ax, sample, (gopts['k_info']['classes']['class_id'] != (i+1)), 'cluster %d' % ((i+1),), color=clust_colors[i%len(clust_colors)])
                
        sys.stderr.write("    -> Generating average profile...\n")
        if gopts['args'].rotate:
            color = gopts['args'].colors[sample.sig_id % len(gopts['args'].colors)]
        else:
            color = gopts['args'].colors[sample.bed_id % len(gopts['args'].colors)]
        
        make_average_sig_plot(ax, sample, color)
        
        if 'kavg' in gopts['args'].plot:
            leg = ax.legend(loc='best', bbox_to_anchor=(0.5, -0.1))
            leg.get_frame().set_linewidth(0.1)
#end add_signal_to_figure()


def make_violin_plot(ax, sample, color='k'):
    
    summary = ttstats.summarize_data(sample.signal_array, method=gopts['args'].summarymethod, axis=1)
    #print summary.shape
    #print gopts['x_axis'].shape
#    label = sample.sig_label if gopts['args'].rotate else sample.bed_label
    if gopts['args'].rotate:
        positions = [sample.sig_id]
    else:
        positions = [sample.bed_id]

    vparts = ax.violinplot(summary, positions=positions, points=len(summary), showmeans=False, showmedians=False, showextrema=False)
    
    #change the filled area of the violin
    for pc in vparts['bodies']:
        pc.set_facecolor(color)
        #pc.set_edgecolor('black')
        pc.set_alpha(1)
    
    quartile1, medians, quartile3 = np.percentile(summary, [25, 50, 75])
    whiskers = adjacent_values(sorted(summary), quartile1, quartile3)
    whiskersMin, whiskersMax = whiskers[0], whiskers[1]
    
    ax.scatter(positions, medians, marker='o', color=color, s=0.2, zorder=3)
    ax.vlines(positions, quartile1, quartile3, color='k', linestyle='-', lw=5, alpha=0.5)
    ax.vlines(positions, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1, alpha=0.5)
    #ax.vlines(positions, medians, medians, color='white', linestyle='-', lw=5)

    ax.set_ylim(bottom=sample.raw_min, top=sample.raw_max)
    ax.set_xticklabels([])
    ax.set_xticks([])
    if not sample.show_yaxis:
        ax.set_yticklabels([])
        
    return ax
#end make_average_sig_plot()

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def make_average_sig_plot(ax, sample, color='k'):
    if gopts['args'].vline:
        ax.axvline(0, linestyle=gopts['args'].vlinestyle, color='k', linewidth=gopts['args'].vlineweight)
        if gopts['args'].align == 'scale':
            ax.axvline(gopts['args'].scaleregionsize, linestyle=gopts['args'].vlinestyle, color='k', linewidth=gopts['args'].vlineweight)
    
    summary = ttstats.summarize_data(sample.signal_array, method=gopts['args'].summarymethod, axis=0)
    #print sample.signal_array.shape
    #print summary.shape
    #print gopts['x_axis'].shape
    label = sample.sig_label if gopts['args'].rotate else sample.bed_label
    ax.plot(gopts['x_axis'], summary, color=color, label=label, alpha=0.7, linewidth=gopts['args'].linewidth)
    
    
    if gopts['args'].showci:
        computed_error = ttstats.compute_error(sample, gopts['args'].ciwidth)
        ax.fill_between(gopts['x_axis'], summary, summary + computed_error, facecolor=color, edgecolor='none', alpha=0.2)
        ax.fill_between(gopts['x_axis'], summary, summary - computed_error, facecolor=color, edgecolor='none', alpha=0.2)
        
    ax.set_xlim(gopts['x_axis'][0], gopts['x_axis'][-1])
    ax.set_ylim(bottom=sample.avg_min, top=sample.avg_max)
    format_tick_marks(ax)
        
    if not sample.show_yaxis:
        ax.set_yticklabels([])
    return ax
#end make_average_sig_plot()


def add_masked_group_to_avg_plot(ax, sample, mask, label, color='k'):
    #print mask
    real_mask = np.zeros(sample.signal_array.shape, dtype=np.bool)
    counts = [0, 0]
    for maskrow in xrange(mask.shape[0]):
        if mask[maskrow]:
            counts[1] = counts[1] + 1
            real_mask[maskrow,] = True
        else:
            counts[0] = counts[0] + 1

    if gopts['args'].showci:
        metaseq.plotutils.ci_plot(gopts['x_axis'], np.ma.array(sample.signal_array, mask=real_mask), gopts['args'].ciwidth, ax, line_kwargs=dict(color=color, label=label), fill_kwargs=dict(color=color, alpha=0.3))
    else:
        ax.plot(gopts['x_axis'], np.ma.array(sample.signal_array, mask=real_mask).mean(axis=0), color=color, label=label)
#end add_masked_group_to_avg_plot

def make_sig_heatmap(ax, sample):
    #sys.stderr.write("-> saturation: (%d, %d)\n" % (vmin, vmax))
    sample.mappable = metaseq.plotutils.imshow(
            sample.signal_array,
            x=gopts['x_axis'],
            ax=ax,
            vmin=sample.heat_min, 
            vmax=sample.heat_max, 
            percentile=False,
            sort_by=sample.sort_order
    )
    format_tick_marks(ax)
    
    #side label options
    if gopts['args'].ylabelrot != 90:
        side_label_opts = {'rotation': gopts['args'].ylabelrot, 'verticalalignment': 'top', 'horizontalalignment': 'right'}
    else:
        side_label_opts = {}
        
    #top label options
    if gopts['args'].xlabelrot != 0:
        top_label_opts = {'rotation': gopts['args'].xlabelrot, 'verticalalignment': 'bottom', 'horizontalalignment': 'left'}
    else:
        top_label_opts = {}
    
    if gopts['args'].rotate:
        if sample.sig_id == 0: #first row
            ax.set_title(sample.bed_label, **top_label_opts)
        if sample.bed_id == 0: #first column
            ax.set_ylabel(sample.sig_label, **side_label_opts)
    else:
        if sample.bed_id == 0: #first row
            ax.set_title(sample.sig_label, **top_label_opts)
        if sample.sig_id == 0: #first column
            ax.set_ylabel(sample.bed_label, **side_label_opts)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if gopts['args'].vline:
        ax.axvline(0, linestyle=gopts['args'].vlinestyle, color='k', linewidth=gopts['args'].vlineweight)
        if gopts['args'].align == 'scale':
            ax.axvline(gopts['args'].scaleregionsize, linestyle=gopts['args'].vlinestyle, color='k', linewidth=gopts['args'].vlineweight)
    if gopts['args'].hline:
        for hl in sample.hlines:
            ax.axhline(hl, linestyle=gopts['args'].hlinestyle, color='k', linewidth=gopts['args'].hlineweight)
#end make_sig_heatmap()

def make_colormap_strip_for_groups(fig, groups, samples):
    groups_list = get_groups(groups, samples)
    for i in range(len(groups_list)):
        make_colormap_strip(fig, [s for s in samples if str(s.sig_id) in groups_list[i]])
#end make_colormap_strip_for_groups()

def make_colormap_strip(fig, samples):
    cmap = colormap_adjust.smart_colormap(samples[0].heat_min, samples[0].heat_max, color_high=gopts['args'].heatcolor)
    norm = mpl.colors.Normalize(vmin=samples[0].heat_min, vmax=samples[0].heat_max)
    cols = []
    for s in samples:
        cols.append(s.sig_id)
    #max_col = max(cols)
    cbar_axis = get_plot_axes('cbar', samples[0].group, min(cols), max(cols))
    cur_pos = cbar_axis.get_position()
    new_pos = [cur_pos.x0, cur_pos.y0, cur_pos.width * 0.25, cur_pos.height]
    cbar_axis.set_position(new_pos)
    [i.set_linewidth(0.1) for i in cbar_axis.spines.itervalues()]
    cb = mpl.colorbar.ColorbarBase(cbar_axis, cmap=cmap, norm=norm)
    cb.outline.set_linewidth(0.1)
#end make_colormap_strip()

def format_tick_marks(ax):
    if gopts['co'].align == 'scale':
        majors = [-gopts['co'].upstream, 0, gopts['co'].scaleregionsize, gopts['co'].scaleregionsize + gopts['co'].downstream]
        labels = []
        for m in majors:
            if m == 0:
                labels.append('S')
            elif m == gopts['co'].scaleregionsize:
                labels.append('E')
            else:
                if m > 0:
                    labels.append('{:+g}'.format((m - gopts['co'].scaleregionsize) / gopts['co'].units[0]))
                else:
                    labels.append('{:+g}'.format(m / gopts['co'].units[0]))
        ax.set_xticks(majors)
        ax.set_xticklabels(labels)
    else:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:+g}'.format(x / gopts['co'].units[0]))) # display with the proper units
#end format_tick_marks()


def save_figure(fig, notes):
    if notes != "":
        notes = "."+notes
    savename = "%s%s.%s" % (gopts['output_base'], notes, gopts['args'].format)
    sys.stderr.write('Saving Figure.....\n')
    fig.savefig(savename, dpi=gopts['args'].dpi, bbox_extra_artists=gopts['extra_artists'], bbox_inches='tight')
    sys.stderr.write(' => See %s for results.\n' % (savename,))
#end save_figure()

if __name__ == "__main__":
    main()
