#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import argparse
import metaseq
import pybedtools
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import multiprocessing
import os
import sys
import collections
from scipy import stats
import fastcluster
import time
import subprocess
from sklearn import metrics
import scipy.cluster.hierarchy as sch
from ThackTech import filetools, chromtools
from ThackTech.Plotting import stats as ttstats, sigcollector





gopts = {
    'args': None,
    'num_bins': 0,
    'output_base': '',
    'fig_cols': 0,
    'fig_rows': 0,
    'x_axis': None,
    'plot_axes': {},
    'extra_artists': []
}


class ProfileSample:

    def __init__(self, id, sig_id, bed_id, signal_array, sig_label, bed_label):
        self.id = id
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
    distance_metrics = ['euclidean','minkowski','cityblock','seuclidean','sqeuclidean','cosine','correlation','hamming','jaccard','chebyshev','canberra','braycurtis','mahalanobis','yule','matching','dice','kulsinski','rogerstanimoto','russellrao','sokalmichener','sokalsneath','wminkowski']
    linkage_methods = ['single','complete','average','weighted','centroid','median','ward']

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bed', action='append', help='Bed files to profile against.')
    parser.add_argument('-s', '--sig', action='append', required='true', help='Signal files to profile. BigWig, BAM, BED, BigBed files are supported, though all but bigwig will be treated by counting intervals.')
    parser.add_argument('-sl', '--slabel', action='append', default=[], help='Signal labels for each plot. Specify in the same order as the signal files. If not supplied, file basename will be used.')
    parser.add_argument('--title', action='store', default='', help='Title for entire plot')
    
    profile_group = parser.add_argument_group('Profiling Options')
    profile_group.add_argument('--scalegroups', action='store', default=None, help='Groups of plots to share color/y-axis scales. If not specified, all plots will be constructed with the same scale. Parameter should be specified as comma-separated lists of 0-based offsets of samples, and groups separated with a semicolon. Ex: 0;1,2;3,4 results in sample 0 plotted with independent scale, 1 and 2 sharing scale, and 3 and 4 sharing scale. If specified, but parameter omits samples, then the omitted samples will each be scaled independently.')
    profile_group.add_argument('--scalemethod', action='store', choices=['max', 'pooled'], default='max', help='Method used when scaling signal groups. "max" will use the maximum of each individual within the group, while "pooled" will pool data first then compute the scale.')
    profile_group.add_argument('--scalemetric', action='store', choices=['max', 'outliermax', 'percentile'], default='outliermax', help='Metric used when scaling signal groups. "max" will use the real maximum, "outliermax" will use the maximum after removing outliers using a z-scoring method (Hoaglin et. al. 1993), and "percentile" will compute the max at the given percentile.')
    default_scalemetricthreshold = {'outliermax': 3.5, 'percentile': 0.01 }
    profile_group.add_argument('--scalemetricthreshold', action='store', default=None, type=float, help='Threshold used for --scalemetric, but ignored when --scalemetric=max. For "outliermax" this is a modified z-score (based on the median absolute deviation) where z-scores greater than this value will be classified as outliers (default %f). For "percentile" this is the top percentile for which values are ignored fo scale computation (default %f)' % (default_scalemetricthreshold['outliermax'], default_scalemetricthreshold['percentile']))
    profile_group.add_argument('--nochromfilter', action='store_true', help='Do not filter bed file for common chromosomes. By default profileing only occurs on chromosomes common to all files (useful for ignoring random or Un* chromosomes).')
    profile_group.add_argument('--up', action='store', default=1000, type=int, help='Span upstream (in bp) from selected intervals alignement to profile.')
    profile_group.add_argument('--down', action='store', default=1000, type=int, help='Span downstream (in bp) from selected intervals alignement to profile.')
    profile_group.add_argument('--align', action='store', choices=['left', 'center', 'right', 'scale'], default='center', help='Method used to align intervals. Left aligns intervals along the 5\' end, while right aligns intervals aling the 3\' end (assuming intervals provide strand information and --dir is specified, otherwise + strand is assumed). Center aligns intervals along the center center point of intervals.')
    profile_group.add_argument('--scaleregionsize', action='store', type=int, default='3000', help='Size of the scaled region (i.e. gene-body) in bp. This option is only useful for "--align scale" option.')
    profile_group.add_argument('--res', action='store', default=10, type=int, help='Profiling resolution in bp.')
    profile_group.add_argument('--dir', action='store_true', help='If set, the direction (+/-) [strand orientation] is considered in profiling. Strand information should be present in bed file used, and if not present is assumed to be +.')
    profile_group.add_argument('--genome', action='store', help='UCSC reference genome that intervals and signal are computed with (i.e. mm9, hg19, etc.).')

    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--name', action='store', default='', help='Base name for the plot output.')
    output_group.add_argument('--format', action='store', default='pdf', choices=['pdf', 'png', 'svg'], help='Format to output the final figure.')
    available_plot_types = ['hist','2dhist','corrmatrix','heat']
    output_group.add_argument('--plot', action='append', choices=available_plot_types+['all'], required='true', help='Types of plots to produce. Supply multiple options to produce hybrid plots. avg will produce average profiles, heat will produce heatmaps, kavg will produce average profiles for each class determined by k-means clustering (only available when used with --plot heat AND --kmeans).')
    output_group.add_argument('--dpi', action='store', type=int, default=300, help='DPI resolution of the saved figure.')
    output_group.add_argument('--width', action='store', type=float, default=8, help='Width (in inches) of the figure.')
    output_group.add_argument('--height', action='store', type=float, default=6, help='Height (in inches) of the figure.')
    output_group.add_argument('--dump', action='store_true', help='Output raw profiles.')
    
    corr_group = parser.add_argument_group('Correlation Plotting Options')
    corr_group.add_argument('--corrshowr', action='store_true', help='For correlation matrix, show the correlation coefficient on the plot itself.')
    corr_group.add_argument('--corrshowp', action='store_true', help='For correlation matrix, show the p-value for the correlation coefficient on the plot itself.')
    corr_group.add_argument('--corrmethod', action='store', choices=['spearman', 'pearson'], default='spearman', help='For correlation matrix, method for computing correlation coefficient.')
    corr_group.add_argument('--corrdendro', action='store_true', help='For correlation matrix, show dendrogram.')
    corr_group.add_argument('--corrlinkagemethod', action='store', choices=linkage_methods, default=linkage_methods[0], help='For dendrogram computation,linkage method to use. See scipy.cluster.hierarchy.linkage documentation for details on methods.')
    corr_group.add_argument('--corrdistancemetric', action='store', choices=distance_metrics, default=distance_metrics[0], help='For dendrogram computation, distance metric to use. See scipy.spatial.distance.pdist documentation for details on metrics.')
    
    heat_group = parser.add_argument_group('Heat Plotting Options')
    heat_group.add_argument('--heatsort', action='store', type=int, default=None, help='0-based index of the sample to sort by for heat plot.')
    heat_group.add_argument('--heatdendro', action='store_true', help='For heat plot, show dendrogram.')
    heat_group.add_argument('--heatlinkagemethod', action='store', choices=linkage_methods, default=linkage_methods[0], help='For dendrogram computation,linkage method to use. See scipy.cluster.hierarchy.linkage documentation for details on methods.')
    heat_group.add_argument('--heatdistancemetric', action='store', choices=distance_metrics, default=distance_metrics[0], help='For dendrogram computation, distance metric to use. See scipy.spatial.distance.pdist documentation for details on metrics.')
    
    hist_group = parser.add_argument_group('Histogram Plotting Options')
    hist_group.add_argument('--histbins', action='store', type=int, default=50, help='Number of bins to use in histogram plots.')
    hist_group.add_argument('--hist2dsmooth', action='store_true', help='Smooth the 2D-Histogram.')
    hist_group.add_argument('--hist2dregression', action='store_true', help='Show a polynomial fit of the data.')
    hist_group.add_argument('--hist2dregdegree', action='store', type=int, default=1, help='Degree of the polynomial to fit the data.')
    hist_group.add_argument('--histlog', action='store_true', help='Log normalize frequency values.')
    
    plot_group = parser.add_argument_group('General Plotting Options')
    plot_group.add_argument('--cmap', action='store', default='YlGnBu', help='Color map to use for plotting.')
    plot_group.add_argument('--fontsize', action='store', type=int, default=8, help='Font size used in plots.')
    plot_group.add_argument('--legendfontsize', action='store', type=int, default=6, help='Font size used in plot legend.')
    plot_group.add_argument('--numticks', action='store', type=int, default=4, help='Number of tick marks to use in the x-axis.')
    
    resources_group = parser.add_argument_group('Resource Options')
    resources_group.add_argument('--cache', action='store_true', help='If set, profiles are dumped to disk.')
    resources_group.add_argument('--cachedir', action='store', default=".correlator_cache", help='Location to place cache files.')
    resources_group.add_argument('-p', '--processors', action='store', dest='cpus', metavar='cpus', default=system_cpu_count, type=int, help='Number of processors to use.')
    gopts['args'] = args = parser.parse_args()
    if 'all' in args.plot:
        args.plot = available_plot_types
    document_args()

    # if len(args.bed) > 1 and len(args.sig) > 1:
        # sys.stderr.write("ERROR: Multiple signal AND bed files are not allowed!\n")
        # sys.exit(1)
    if len(args.slabel) == 0:
        sys.stderr.write("WARNING: Signal labels were not supplied. Using signal file basenames instead.\n")
    # if len(args.ilabel) == 0:
        # sys.stderr.write("WARNING: Interval labels were not supplied. Using interval file basenames instead.\n")
    if len(args.slabel) > 0 and len(args.slabel) < len(args.sig):
        sys.stderr.write("ERROR: Fewer signal labels supplied than there are signals! Signal and label counts must be equal if supplying labels!\n")
        sys.exit(1)
    # if len(args.ilabel) > 0 and len(args.ilabel) < len(args.bed):
        # sys.stderr.write("ERROR: Fewer interval labels supplied than there are intervals! Interval and label counts must be equal if supplying labels!\n")
        # sys.exit(1)
    # if len(args.inp) > 0 and len(args.inp) < len(args.sig):
        # sys.stderr.write("ERROR: Fewer inputs supplied than there are samples! sample and input counts must be equal if supplying inputs!\n")
        # sys.exit(1)
    if len(args.plot) <= 0:
        sys.stderr.write("ERROR: No plot types were selected!\n")
        sys.exit(1)
    if args.cpus <= 0:
        sys.stderr.write("WARNING: Cannot use fewer than one processor for computation. Setting --processors = 1!\n")
        args.cpus = 1
    if args.cpus > system_cpu_count:
        sys.stderr.write("WARNING: Cannot use more than the number of processors on the system ("+system_cpu_count+") for computation. Setting --processors = "+system_cpu_count+"!\n")
        args.cpus = system_cpu_count
    # if len(args.bed) > 1 and args.kmeans:
        # sys.stderr.write("ERROR: Cannot perform K-means analysis when using multiple beds!\n")
        # sys.exit(1)
    # if 'kavg' in args.plot and not args.kmeans and not 'heat' in args.plot:
        # sys.stderr.write("WARNING: to plot averages for clusters, --kmeans and --plot heat must be used! Ignoring this directive!\n")
        # args.plot.remove('kavg')
    
    
    
    collection_opts = sigcollector.CollectorOptions()
    collection_opts.align = args.align
    collection_opts.upstream = abs(args.up)
    collection_opts.downstream = abs(args.down)
    collection_opts.scaleregionsize = abs(args.scaleregionsize)
    collection_opts.resolution = abs(args.res)
    collection_opts.direction = bool(args.dir)
    collection_opts.validate()
    
    
    
    gopts['output_base'] = "%s.%du_%dd_%s" % (args.name, args.up, args.down, (args.align if args.align == 'scale' else args.align))
    gopts['savename_notes'] = []
    
    
    
    
    chromsets = chromtools.get_common_chroms(args.bed, args.sig)
    if not args.nochromfilter:
        sys.stderr.write("Filtering common chromomsomes....\n")
        chromsets.use = chromsets.common
        sys.stderr.write("=> Keeping:  "+', '.join(map(str,sorted(chromsets.common)))+"\n")
        sys.stderr.write("=> Ignoring: "+', '.join(map(str,sorted(chromsets.uncommon)))+"\n")
        sys.stderr.write("\n")
    else:
        chromsets.use = chromsets.all
    
    samples = []
    for s in xrange(len(args.sig)):
        for b in xrange(len(args.bed)):
            sys.stderr.write("Processing %s vs %s....\n" % (args.bed[b], args.sig[s]))
            sys.stderr.write("-> Preparing intervals.....\n")
            
            bedtool = sigcollector.IntervalProvider(args.bed[b], collection_opts, args.genome, gopts['chromsets'].use)
            s_label = args.slabel[s] if len(args.slabel)-1 >= s else os.path.splitext(os.path.basename(args.sig[s]))[0]
            b_label = os.path.splitext(os.path.basename(args.bed[b]))[0]#args.ilabel[b] if len(args.ilabel)-1 >= b else os.path.splitext(os.path.basename(args.bed[b]))[0]

            signal = sigcollector.get_signal(bedtool, s_label+b_label, args.sig[s], None, cache_dir=(args.cachedir if args.cache else None), cache_base=args.name, collectionmethod=args.collectionmethod, cpus=args.cpus)
            samples.append(ProfileSample(len(samples), s, b, signal, s_label, b_label))
            sys.stderr.write("\n")
    
    gopts['group_count'] = 0 #count_groups(args.scalegroups, samples)
    gopts['fig_cols'] = len(args.sig) + gopts['group_count']
    gopts['fig_rows'] = len(args.sig) #((args.heatplotrows if 'heat' in args.plot else 0) * len(args.bed)) + (args.avgplotrows if 'avg' in args.plot else 0)
    plt.rcParams['font.size'] = args.fontsize
    plt.rcParams['legend.fontsize'] = args.legendfontsize

    #compute saturation points and average profile limits
    if (args.scalemetric in default_scalemetricthreshold) and (args.scalemetricthreshold is None):
        args.scalemetricthreshold = default_scalemetricthreshold[args.scalemetric]
    compute_group_scales(args.scalegroups, samples)
    
    
    
    if 'hist' in args.plot:
        make_1d_hist_plots(samples)
    if '2dhist' in args.plot:
        make_2d_hist_plots(samples)
    if 'corrmatrix' in args.plot:
        make_correlation_plot(samples)
    if 'heat' in args.plot:
        make_heat_plot(samples)
    
    if args.dump:
        for b in xrange(len(args.bed)):
            bedtool = sigcollector.IntervalProvider(args.bed[b], collection_opts, args.genome, gopts['chromsets'].use)
            dump_raw_data([s for s in samples if (s.bed_id) == b], bedtool)
            
    sys.stderr.write('Done!')
#end main()

def document_args():
    sys.stderr.write("\n")
    sys.stderr.write("# ARGUMENTS:\n")
    sys.stderr.write("# ==================================================\n")
    
    sys.stderr.write("# Run Start: %s\n" % (time.strftime("%a, %d %b %Y %H:%M:%S"),))
    sys.stderr.write("# Name: %s\n" % (gopts['args'].name,))
    sys.stderr.write("# Interval File(s):\n")
    for i in range(len(gopts['args'].bed)):
        sys.stderr.write("#\t%d: %s\n" % (i, os.path.abspath(gopts['args'].bed[i]),))
    sys.stderr.write("# Signal File(s):\n")
    for i in range(len(gopts['args'].sig)):
        sys.stderr.write("#\t%d: %s\n" % (i, os.path.abspath(gopts['args'].sig[i]),))
        
    sys.stderr.write("# --------------------------------------------------\n")
    sys.stderr.write("# Alignment: %s\n" % (gopts['args'].align,))
    if gopts['args'].align == 'scale':
        sys.stderr.write("# Span: %dbp Up; %dbp Down from scaled interval\n" % (gopts['args'].up, gopts['args'].down))
    else:
        sys.stderr.write("# Span: %dbp Up; %dbp Down\n" % (gopts['args'].up, gopts['args'].down))
    sys.stderr.write("# Filter common chromosomes: %s\n" % ('OFF' if gopts['args'].nochromfilter else 'ON',))
    sys.stderr.write("# Scale Groups: %s\n" % (gopts['args'].scalegroups,))
    #groups = get_groups(gopts['args'].scalegroups
    #for g i
    
    
    #sys.stderr.write("# Saturate Min: %f\n" % (gopts['args'].saturatemin,))
    #sys.stderr.write("# Saturate Max: %f\n" % (gopts['args'].saturatemax,))
    sys.stderr.write("# Plot Types: %s\n" % (", ".join(gopts['args'].plot),))
    sys.stderr.write("# Cache Data: %s\n" % ('ON' if gopts['args'].cache else 'OFF',))
    if gopts['args'].cache:
        sys.stderr.write("# Cache Location: %s\n" % (os.path.abspath(gopts['args'].cachedir),))
    sys.stderr.write("# Num Threads: %d\n" % (gopts['args'].cpus,))
    if 'hist' in gopts['args'].plot:
        sys.stderr.write("# --------------------------------------------------\n")
        sys.stderr.write("# 1D-Histogram Options:\n")
        sys.stderr.write("# Number of Bins: %d\n" % (gopts['args'].histbins,))
        sys.stderr.write("# Log Scale Histogram: %s\n" % ('OFF' if gopts['args'].histlog else 'ON',))
        
    if '2dhist' in gopts['args'].plot:
        sys.stderr.write("# --------------------------------------------------\n")
        sys.stderr.write("# 2D-Histogram Options:\n")
        sys.stderr.write("# Number of Bins: %d\n" % (gopts['args'].histbins,))
        sys.stderr.write("# Log Scale Histogram: %s\n" % ('OFF' if gopts['args'].histlog else 'ON',))
        
        sys.stderr.write("# Smooth 2D-Histogram: %s\n" % ('OFF' if gopts['args'].hist2dsmooth else 'ON',))
        sys.stderr.write("# Show a Polynomial Fit: %s\n" % ('OFF' if gopts['args'].hist2dregression else 'ON',))
        if gopts['args'].hist2dregression:
            sys.stderr.write("# Degree of polynomial fit: %d\n" % (gopts['args'].hist2dregdegree,))
            
    if 'corrmatrix' in gopts['args'].plot:
        sys.stderr.write("# --------------------------------------------------\n")
        sys.stderr.write("# Correlation Matrix Options:\n")
        sys.stderr.write("# Correlation Method: %s\n" % (gopts['args'].corrmethod,))
        sys.stderr.write("# Show Correlation Coefficient on Plot: %s\n" % ('OFF' if gopts['args'].corrshowr else 'ON',))
        sys.stderr.write("# Show correlation p-value on Plot: %s\n" % ('OFF' if gopts['args'].corrshowp else 'ON',))
        sys.stderr.write("# Show Dendrogram: %s\n" % ('OFF' if gopts['args'].corrdendro else 'ON',))
        if gopts['args'].corrdendro:
            sys.stderr.write("# Linkage Method: %s\n" % (gopts['args'].corrlinkagemethod,))
            sys.stderr.write("# Distance Metric: %s\n" % (gopts['args'].corrdistancemetric,))
        
    if 'heat' in gopts['args'].plot:
        pass
        #make_heat_plot(samples)
    
    
    
    sys.stderr.write("\n")
#end document_args()

def get_groups(groups, samples):
    if groups is None:
        return [list(set([str(s.sig_id) for s in samples]))]
    else:
        non_covered_samples = set([str(s.sig_id) for s in samples]) - set(groups.replace(";",",").split(","))
        if len(non_covered_samples) > 0:
            groups += ";" + ";".join(non_covered_samples)
        final_groups = []
        for group in groups.split(";"):
            final_groups.append([str(s.sig_id) for s in samples if str(s.sig_id) in group.split(",")])
        return final_groups
#end get_groups()

def count_groups(groups, samples):
    return len(get_groups(groups, samples))
#end count_groups()

def compute_group_scales(groups, samples):
    groups_list = get_groups(groups, samples)
    for i in range(len(groups_list)):
        compute_scales_for_group(i, [s for s in samples if str(s.sig_id) in groups_list[i]])
#end compute_group_scales()

def compute_scales_for_group(group_id, samples):
    if gopts['args'].scalemethod == 'pooled':
        if len(samples) > 1:
            pooled = np.vstack(tuple([s.signal_array for s in samples]))
        else:
            pooled = samples[0].signal_array
        pooled_min, pooled_max = compute_data_scale(pooled)
    else:
        pooled_max = sys.float_info.min
        pooled_min = 0
        for s in samples:
            localmin, localmax = compute_data_scale(s.signal_array)
            pooled_min = min(pooled_min, localmin)
            pooled_max = max(pooled_max, localmax)
    
    pooled_min = 0
    for s in samples:
        s.group = group_id
        s.min = pooled_min
        s.max = pooled_max
    sys.stderr.write("-> Computed min/max for group %d (%f, %f)\n" % (group_id, pooled_min, pooled_max))
#end compute_scales_for_group()


def compute_data_scale(data):
    if gopts['args'].scalemetric == 'max':
        return (np.min(data), np.max(data))
    elif gopts['args'].scalemetric == 'percentile':
        return (np.percentile(data, (gopts['args'].scalemetricthreshold*100)), np.percentile(data, ((1-gopts['args'].scalemetricthreshold)*100)))
    else: #outliermax
        use = data[(~ttstats.is_outlier(data, gopts['args'].scalemetricthreshold))]
        return (np.min(use), np.max(use))
#end compute_data_scale()


# def get_signal(bedtool, sig_file, label):
#     cache_dir = os.path.abspath(gopts['args'].cachedir)
#     if not os.path.exists(cache_dir):
#         os.makedirs(cache_dir)
#     cache_name = os.path.join(cache_dir, "%s.%s" % (gopts['output_base'], filetools.make_label_filename_safe(label)))
#     
#     if not gopts['args'].cache or not os.path.exists(cache_name + '.npz'):
#         sys.stderr.write("-> Loading signal....\n")
#         sig = metaseq.genomic_signal(sig_file, detect_signal_type(sig_file))
#         sys.stderr.write("-> Computing signal at intervals....\n")
#         sig_array = sig.array(bedtool, bins=1, method='summarize', processes=gopts['args'].cpus)
#             
#         if gopts['args'].cache:
#             sys.stderr.write("-> Persisting data to disk...\n")
#             cache_data = {label: sig_array}
#             metaseq.persistence.save_features_and_arrays(features=bedtool,
#                                                          arrays=cache_data,
#                                                          prefix=cache_name,
#                                                          #link_features=True,
#                                                          overwrite=True)
#     else:
#         sys.stderr.write("-> Loding data from cache....\n")
#         features, arrays = metaseq.persistence.load_features_and_arrays(prefix=cache_name)
#         sig_array = arrays[label]
#     
#     return sig_array
# #end get_signal()


def dump_raw_data(samples, bed):
    intervals = list(bed)
    with open(generate_save_name('tsv', 'dump'), "w") as file:
        file.write('Chr\tStart\tStop\tName\tBedScore\tStrand')
        for s in samples:
            file.write('\t%s' % (s.sig_label,))
        if 'heat' in gopts['args'].plot and gopts['args'].heatsort is not None:
            file.write('\tHeatOrder')
        file.write('\n')#end header row
        
        for i in range(len(samples[0].signal_array)):
            file.write('%s\t%d\t%d\t%s\t%s\t%s' % (intervals[i].chrom, intervals[i].start, intervals[i].stop, intervals[i].name, str(intervals[i].score), intervals[i].strand))
            for s in samples:
                file.write('\t%f' % (s.signal_array[i],))
            if 'heat' in gopts['args'].plot and gopts['args'].heatsort is not None:
                file.write('\t%d' % (samples[0].heatsort[i],))
            file.write('\n')
#end dump_raw_data()


def compute_sorting(samples, sort_index, method):
    if sort_index is None:
        for s in samples:
            s.sort_order = None
    else:
        sort_orders = {} #sort orders indexed by interval id
        for s in samples:
            if s.sig_id == sort_index:
                sort_orders[s.bed_id] = getattr(s.signal_array, method)(axis=1)
        for s in samples:
            s.sort_order = sort_orders[s.bed_id]
#end compute_sorting()

class PlotPosition:
    def __init__(self, row, col, row_span, col_span):
        self.row = row
        self.col = col
        self.rs = row_span
        self.cs = col_span
        self.pos = (row, col)
    #end __init__
#end class PlotPosition


def make_heat_plot(samples):
    
    
    num_samples = len(samples)
    num_rows = num_samples + 1
    num_cols = num_samples + 2
    heat_loc = PlotPosition(1, 1, num_samples, num_samples)
    tden_loc = PlotPosition(0, 1, 1, num_samples)
    lden_loc = PlotPosition(1, 0, num_samples, 1)
    cbar_loc = PlotPosition(1, num_cols-1, num_samples, 1)
    
    data = np.hstack(tuple([(s.signal_array / s.max) for s in samples]))
    
    fig = plt.figure(figsize=(9,8), dpi=gopts['args'].dpi)
    
    if gopts['args'].heatsort is not None:
        order = np.argsort(data[:,gopts['args'].heatsort])
        data = data[order,:]
        samples[0].heatsort = np.argsort(order)
    
    if gopts['args'].heatdendro:
        #top dendrogram
        ax1 = plt.subplot2grid((num_rows,num_cols), tden_loc.pos, rowspan=tden_loc.rs, colspan=tden_loc.cs)
        Y = fastcluster.linkage(np.rot90(data), method=gopts['args'].heatlinkagemethod, metric=gopts['args'].heatdistancemetric)
        data = data[:,sch.leaves_list(Y)]#reorder the correlation matrix to reflect dendrogram ordering
        Z1 = sch.dendrogram(Y, labels=[s.sig_label for s in samples], leaf_font_size=8)
        ax1.set_yticks([])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        
        if gopts['args'].heatsort is None:
            #left dendrogram
            ax2 = plt.subplot2grid((num_rows,num_cols), lden_loc.pos, rowspan=lden_loc.rs, colspan=lden_loc.cs)
            Y = fastcluster.linkage(data, method=gopts['args'].heatlinkagemethod, metric=gopts['args'].heatdistancemetric)
            data = data[sch.leaves_list(Y),:]#reorder the correlation matrix to reflect dendrogram ordering
            sch.dendrogram(Y, orientation='right', p=5, truncate_mode='lastp', leaf_rotation=90, leaf_font_size=8)
            ax2.set_xticks([])
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.spines['left'].set_visible(False)
    
    axheat = plt.subplot2grid((num_rows,num_cols), heat_loc.pos, rowspan=heat_loc.rs, colspan=heat_loc.cs)
    heatmap = axheat.pcolormesh(data, vmax=1, cmap=plt.cm.get_cmap(gopts['args'].cmap))
    axheat.set_frame_on(False)
    axheat.xaxis.tick_top()
    axheat.set_xticklabels([s.sig_label for s in samples], minor=False)
    axheat.set_yticks([])
    axheat.set_xticks([])
    axheat.set_ylim(0, data.shape[0])
    
    # Plot colorbar.
    axcolor = plt.subplot2grid((num_rows,num_cols), cbar_loc.pos, rowspan=cbar_loc.rs, colspan=cbar_loc.cs)
    plt.colorbar(heatmap, cax=axcolor)
    axcolor.set_ylabel("Relative Signal")

    save_close_figure(fig, ['heat', gopts['args'].heatlinkagemethod, gopts['args'].heatdistancemetric])
#end make_heat_plot()




def make_correlation_plot(samples):

    corr_matrix = compute_correlation_matrix(samples, gopts['args'].corrmethod)
    
    num_samples = len(samples)
    num_rows = len(samples) + 2
    num_cols = len(samples) + 3
    labels = [s.sig_label for s in samples]

    if gopts['args'].corrdendro:
        # Compute and plot first dendrogram.
        fig = plt.figure(figsize=(9,8), dpi=gopts['args'].dpi)
        ax1 = plt.subplot2grid((num_rows,num_cols), (2, 0), rowspan=num_samples, colspan=2)
        Y = sch.linkage(corr_matrix[:,:,0], method=gopts['args'].corrlinkagemethod, metric=gopts['args'].corrdistancemetric)
        Z1 = sch.dendrogram(Y, orientation='right', labels=labels, leaf_font_size=8)
        ax1.set_xticks([])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        # Compute and plot second dendrogram.
        ax2 = plt.subplot2grid((num_rows,num_cols), (0, 2), rowspan=2, colspan=num_samples)
        Y = sch.linkage(corr_matrix[:,:,0], method=gopts['args'].corrlinkagemethod, metric=gopts['args'].corrdistancemetric)
        Z2 = sch.dendrogram(Y, labels=labels, leaf_rotation=90, leaf_font_size=8)
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        
        #reorder the correlation matrix to reflect dendrogram ordering
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        corr_matrix = corr_matrix[idx1,:]
        corr_matrix = corr_matrix[:,idx2]

    # Plot distance matrix.
    axmatrix = plt.subplot2grid((num_rows,num_cols), (2, 2), rowspan=num_samples, colspan=num_samples)
    im = axmatrix.matshow(corr_matrix[:,:,0], aspect='auto', origin='lower', cmap=plt.cm.get_cmap(gopts['args'].cmap), interpolation='nearest')
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = plt.subplot2grid((num_rows,num_cols), (2, num_samples+2), rowspan=num_samples, colspan=1)
    plt.colorbar(im, cax=axcolor)
    axcolor.set_ylabel("%s Correlation" % (gopts['args'].corrmethod.capitalize(),))
    
    #annotate values
    if gopts['args'].corrshowr or gopts['args'].corrshowp:
        for i in range(corr_matrix.shape[0]):
            for j in range(corr_matrix.shape[1]):
                z = corr_matrix[i,j,0]
                color = im.cmap(im.norm(z))
                if all(c > 0.5 for c in color):
                    color = (0.0, 0.0, 0.0)
                else:
                    color = (1.0, 1.0, 1.0)
                if gopts['args'].corrshowr and gopts['args'].corrshowp:
                    text = '{:0.3f}\n({:.2E})'.format(z, corr_matrix[i,j,1])
                elif gopts['args'].corrshowp:
                    text = '({:.2E})'.format(corr_matrix[i,j,1])
                else:
                    text = '{:0.3f}'.format(z)
                axmatrix.text(j, i, text, ha='center', va='center', color=color)
    
    plt.tight_layout()
    save_close_figure(fig, ['corrmatrix', gopts['args'].corrmethod, gopts['args'].corrlinkagemethod, gopts['args'].corrdistancemetric])
#end make_correlation_plot()


def compute_correlation_matrix(samples, method):
    num_samples = len(samples)
    matrix = np.zeros((num_samples, num_samples), dtype=(float,2))
    
    with open(generate_save_name('tsv', ['corr', method]), "w") as file:
        file.write('sample_1\tsample_2\t%s R\tp\n' % (method,))
        for i in range(num_samples):
            for j in range(num_samples):
                if method == 'spearman':
                    matrix[i,j] = stats.spearmanr(samples[i].signal_array, samples[j].signal_array, axis=None)
                else:
                    matrix[i,j] = stats.pearsonr(samples[i].signal_array, samples[j].signal_array)
                file.write("%s\t%s\t%f\t%f\n" % (samples[i].sig_label, samples[j].sig_label, matrix[i,j][0], matrix[i,j][1]))
    return matrix
#end compute_correlation_matrix()


def make_1d_hist_plots(samples):
    fig = plt.figure(figsize=(8,8), dpi=gopts['args'].dpi)
    num_samples = len(samples)
    num_rows = num_samples
    num_cols = 1
    
    if gopts['args'].histlog:
        norm = mpl.colors.LogNorm()
    else:
        norm = mpl.colors.Normalize(vmin=0.0,vmax=1.0)
    
    for i in range(num_samples):#row
        ax = plt.subplot2grid((num_rows,num_cols), (i, 0), colspan=num_cols)
        normsample1 = samples[i].signal_array.ravel()
        ax.hist(normsample1, gopts['args'].histbins, normed=norm, range=[samples[i].min, samples[i].max])
        ax.set_ylabel(samples[i].sig_label)
    
    save_close_figure(fig, ['hist', '%d-bins' % (gopts['args'].histbins,)])
#end make_1d_hist_plots()


def make_2d_hist_plots(samples):
    
    
    fig_size = (8,9)
    fig = plt.figure(figsize=fig_size, dpi=gopts['args'].dpi)
    num_samples = len(samples)
    num_rows = len(samples) *2  #(exclude diagional)
    num_cols = len(samples) * 2 - 1 #(exclude diagional, but include colorbar)
    smooth = True
    notes = ['2dhist', '%d-bins' % (gopts['args'].histbins,)]
    if gopts['args'].hist2dregression:
        notes.append('%d-orderfit' % (gopts['args'].hist2dregdegree,))
    
    histograms = {}
    
    if gopts['args'].histlog:
        notes.append('lognorm')
        norm = mpl.colors.LogNorm()
    else:
        norm = mpl.colors.Normalize(vmin=0.0,vmax=1.0)
    
    with open(generate_save_name('tsv', notes), "w") as file:
        file.write('sample_1\tsample_2\tmodel\tr^2\n')
        for i in range(num_samples):#row
            for j in range(i):#col
                normsample_y = samples[i].signal_array.ravel() #to be plotted on the y-axis
                normsample_x = samples[j].signal_array.ravel() #to be plotted on the x-axis
                #print samples[i].sig_label, samples[i].signal_array.mean(axis=0)
                
                #compute the polyfit
                coefficients = np.polyfit(normsample_x, normsample_y, gopts['args'].hist2dregdegree)
                polynomial_fit = np.poly1d(coefficients)
                r_squared = metrics.r2_score(normsample_y, polynomial_fit(normsample_x))
                
                file.write("%s\t%s\t%s\t%f\n" % (samples[j].sig_label, samples[i].sig_label, ttstats.format_poly_equation(coefficients), r_squared))
                
                #plot the 2D-histogram with colorbar
                ax = plt.subplot2grid((num_rows,num_cols), ((i-1)*2, j*2), rowspan=2, colspan=2)
                H, xedges, yedges, img = ax.hist2d(normsample_x, normsample_y, bins=gopts['args'].histbins, range=[[samples[j].min, samples[j].max],[samples[i].min, samples[i].max]], norm=norm, cmap=plt.cm.get_cmap(gopts['args'].cmap))#, alpha=0.1, s=10, linewidths=0)
                #plot the polyfit
                xlin = np.linspace(0, samples[j].max, gopts['args'].histbins)
                if gopts['args'].hist2dregression:
                    ax.plot(xlin, polynomial_fit(xlin), '-', c='k')
                
                histograms[(i-1, j)] = {
                    'H': np.rot90(H),
                    'fit_x': xlin,
                    'fit_y': polynomial_fit(xlin)
                }
                
                #plt.colorbar(H[3], ax=ax)
                #make it pretty and informative
                #ax.set_title(sample1.sig_label + ' vs. ' + sample2.sig_label)
                if i == num_samples-1:
                    ax.set_xlabel(samples[j].sig_label)
                    for label in ax.get_xticklabels(): 
                        label.set_rotation(90) 
                else:
                    ax.set_xticklabels([])
                    
                if j == 0:
                    ax.set_ylabel(samples[i].sig_label)
                else:
                    ax.set_yticklabels([])
        
        # Plot colorbar.
        axcolor = plt.subplot2grid((num_rows,num_cols), (0, num_samples-1), rowspan=(num_samples-1)*2, colspan=1)
        plt.colorbar(img, cax=axcolor, norm=norm)
        axcolor.set_ylabel("Relative Frequency")
        
        
        
    if gopts['args'].hist2dsmooth:
        notes.append('smooth')
        plt.clf()
        plt.close()
        fig = plt.figure(figsize=fig_size, dpi=gopts['args'].dpi)
        
        for i in range(num_samples):#row
            for j in range(i):#col
                item = histograms[(i-1, j)]
                ax = plt.subplot2grid((num_rows,num_cols), ((i-1)*2, j*2), rowspan=2, colspan=2)
                img = ax.imshow(item['H'], origin="upper", interpolation="gaussian", aspect='auto', extent=[samples[j].min, samples[j].max, samples[i].min, samples[i].max], norm=norm, cmap=plt.cm.get_cmap(gopts['args'].cmap))
                ax.set_ylim(samples[i].min, samples[i].max)
                ax.set_xlim(samples[j].min, samples[j].max)
                if gopts['args'].hist2dregression:
                    ax.plot(item['fit_x'], item['fit_y'], '-', c='k')
                
                #plt.colorbar(H[3], ax=ax)
                #make it pretty and informative
                #ax.set_title(sample1.sig_label + ' vs. ' + sample2.sig_label)
                if i == num_samples-1:
                    ax.set_xlabel(samples[j].sig_label)
                    for label in ax.get_xticklabels(): 
                        label.set_rotation(90) 
                else:
                    ax.set_xticklabels([])
                    
                if j == 0:
                    ax.set_ylabel(samples[i].sig_label)
                else:
                    ax.set_yticklabels([])
        # Plot colorbar.
        axcolor = plt.subplot2grid((num_rows,num_cols), (0, num_cols-1), rowspan=(num_samples-1)*2, colspan=1)
        plt.colorbar(img, cax=axcolor, norm=norm)
        axcolor.set_ylabel("Relative Frequency")

    
    save_close_figure(fig, notes)
#end make_2d_hist_plots()



def generate_save_name(format, notes=None):
    if notes is not None and hasattr(notes, '__iter__'):
        notes = "_".join(notes)
    if notes != "":
        notes = "."+notes
    savename = "%s%s.%s" % (gopts['output_base'], notes, format)
    return savename
#end generate_save_name()

def save_close_figure(fig, notes=None):
    savename = generate_save_name(gopts['args'].format, notes)
    sys.stderr.write('Saving Figure.....\n')
    fig.savefig(savename, dpi=gopts['args'].dpi, bbox_extra_artists=gopts['extra_artists'], bbox_inches='tight')
    gopts['extra_artists'] = []
    plt.close(fig)
    sys.stderr.write(' => See %s for results.\n' % (savename,))
#end save_figure()

if __name__ == "__main__":
    main()
