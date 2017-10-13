#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from scipy import stats
import argparse
import pandas as pd
from ThackTech.Plotting import stats as ttstats
from ThackTech.Plotting import intensity_plots as ip
from ThackTech import filetools
from ThackTech.Common import make_safe_worker
import itertools
from multiprocessing import Pool

mpl.rcParams['agg.path.chunksize'] = 10000
mpl.rcParams['pdf.fonttype'] = 42 # use TrueType fonts


def main():
	parent_parser = argparse.ArgumentParser(add_help=False)  
	
	stat_choices = ['mean', 'median', 'count', 'sum', 'min', 'max']
	parent_parser.add_argument('--stat', default=stat_choices[0], choices=stat_choices, help="The statistic to compute.")
	parent_parser.add_argument('--out', '-o', action='store', help="Output basename") 
	parent_parser.add_argument('--figformat', '-ff', action='store', choices=['pdf', 'png', 'svg', 'ps', 'eps'], default='pdf', help="Output figure format.") 
	parent_parser.add_argument('--processors', '-p', type=int, default=1, action='store', help="Number of processors to utilize")
	parent_parser.add_argument('--colors', action='store', nargs='+', default=ip.color_cycle, help="Colors to cycle through for plots.")
	parser = argparse.ArgumentParser(add_help=False) 
	subparsers = parser.add_subparsers()

	# create the parser for the "preprocess" command
	preprocess_parser = subparsers.add_parser('preprocess', help='Step #1: Preprocess raw data.', parents=[parent_parser])
	preprocess_parser.add_argument('files', nargs='+', help="Raw data files to process with column 1 as radial, column 2 as intensity")
	preprocess_parser.add_argument('--min', type=float, default=0, help="Lower range of bins.")
	preprocess_parser.add_argument('--max', type=float, default=1, help="Upper range of bins.")
	preprocess_parser.add_argument('--bins', type=int, default=20, help="The number of equal-width bins in the given range.")
	preprocess_parser.add_argument('--rmoutliers', action='store_true', help="When normalizing data, first remove outliers from the calculation of the correction factor.")
	preprocess_parser.add_argument('--skipnorm', action="store_true", help="Do not normalize the data.")
	preprocess_parser.add_argument('--dir', action="store", help="Directory to place results.")
	preprocess_parser.set_defaults(func=preprocess_files)
	
	# create the parser for the "compare" command
	compare_parser = subparsers.add_parser('compare', help='Step #2: Compare pre-processed binned data.', parents=[parent_parser])
	#stat_choices = ['mean', 'median', 'count', 'sum', 'min', 'max']
	#preprocess_parser.add_argument('--stat', default=stat_choices[0], choices=stat_choices, help="The statistic to compute.")
	compare_parser.add_argument('files', nargs='+', help="Raw data files to process with column 1 as radial, column 2 as intensity")
	compare_parser.set_defaults(func=compare_binned_data)
	
	compraw_parser = subparsers.add_parser('compraw', help='Step #2: Compare pre-processed raw data.', parents=[parent_parser])
	#stat_choices = ['mean', 'median', 'count', 'sum', 'min', 'max']
	#preprocess_parser.add_argument('--stat', default=stat_choices[0], choices=stat_choices, help="The statistic to compute.")
	compraw_parser.add_argument('--filegroup', nargs='+', action='append', help="Raw data files to process with column 1 as radial, column 2 as intensity")
	compraw_parser.add_argument('--groupnames', nargs='+', help="Group Names")
	compraw_parser.set_defaults(func=compare_raw_data)
	
	compraw_parser = subparsers.add_parser('coloc', help='Step #2: Compare raw data for colocalization.', parents=[parent_parser])
	#stat_choices = ['mean', 'median', 'count', 'sum', 'min', 'max']
	#preprocess_parser.add_argument('--stat', default=stat_choices[0], choices=stat_choices, help="The statistic to compute.")
	compraw_parser.add_argument('--filegroup', nargs='+', action='append', help="Raw data files to process with column 1 as radial, columns 2 and 3 as intensity")
	compraw_parser.add_argument('--groupnames', nargs='+', help="Group Names")
	compraw_parser.set_defaults(func=coloc_raw_data)
	
	args = parser.parse_args()
	ip.color_cycle = args.colors
	try:
		args.func(args)
	except KeyboardInterrupt:
		sys.stderr.write('Got Keyboard Interrupt.\nGoodbye.\n')
		return
	
#end main()

def preprocess_files(args):
	filetools.ensure_dir(args.dir)
	#worker pool preprocess each data file in parallel
	#Worker performs normalization (if requested) and plotting of raw data
	#worker returns a tuple of filename and dataframe with binned data
	worker_pool = Pool(processes=args.processors)
	results = []
	try:
		for f in args.files:
			worker_pool.apply_async(preprocess_worker, (f, args), callback=results.append)
		worker_pool.close()
		worker_pool.join()
	except KeyboardInterrupt:
		worker_pool.terminate()
		raise
   
	#for f in args.files:
	#	results.append(preprocess_worker(f, args))
	
	#pool together all the results into a master data frame
	master_dfs = {} 
	for i in range(len(results)):
		f = results[i][0]
		data = results[i][1]
		for j in range(1, len(data.columns)):
			if data.columns[j] not in master_dfs:
				master_dfs[data.columns[j]] = pd.DataFrame()
				master_dfs[data.columns[j]]['Radial'] = data['Radial']
			master_dfs[data.columns[j]][get_label_from_path(f)] = data[data.columns[j]]
	
	#write out the pooled data
	for k in master_dfs.keys():
		master_dfs[k].to_csv(os.path.join(args.dir, args.out+'.'+k+'.tsv'), sep='\t', index=False)
#end preprocess_files()

def preprocess_worker(f, args):
	try:
		sys.stderr.write('Processing file {}....\n'.format(f))
		dest_base = os.path.join(args.dir, get_label_from_path(f))
		data = pd.read_csv(f, sep='\t', comment='#', skip_blank_lines=True)
		data = data.dropna()
		if not args.skipnorm:
			dest_base += '.norm'
			data = normalize_data(data, [data.columns[i] for i in range(1, len(data.columns))], args.rmoutliers)
			data.to_csv(dest_base+'.tsv', sep='\t', index=False)
		#print data
		#print data.describe()
		fig = plot_raw_data([data])
		save_figure(fig, dest_base, args.figformat)
			
		#save binned data
		binned_df = bin_dataframe(data, args)
		binned_df.to_csv(dest_base+'.binned.tsv', sep='\t', index=False)
		return (f, binned_df)
	except KeyboardInterrupt:
		return
#end preprocess_worker()

def bin_dataframe(df, args):
	binned_df = pd.DataFrame()
	binned_df['Radial'] = np.linspace(args.min, args.max, args.bins+1, endpoint=True)[1:]
	for i in range(1, len(df.columns)):
		data_stats, bin_edges, binnumber = stats.binned_statistic(df[df.columns[0]], df[df.columns[i]], statistic=args.stat, bins=args.bins, range=[(args.min, args.max)])
		binned_df[df.columns[i]] = data_stats
	return binned_df
#end bin_dataframe()


def fetch_raw_data(args):
	dfs = []
	for i in range(len(args.groupnames)):
		gname = args.groupnames[i]
		cache_name = args.out+'.'+filetools.make_str_filename_safe(gname)+'.pkl.gz'
		if os.path.exists(cache_name):
			print "Loading data for group {}".format(gname)
			dfs.append(pd.read_pickle(cache_name))
		else:
			print "Reading data for group {}".format(gname)
			df = pd.concat([pd.read_csv(f, sep='\t', comment='#', skip_blank_lines=True) for f in args.filegroup[i]])
			df.reset_index(inplace=True, drop=True)
			dfs.append(df)
			df.to_pickle(cache_name)
	return dfs
#end fetch_raw_data()


def compare_raw_data(args):
	dfs = fetch_raw_data(args)
	
	# binned_dfs = []
	# for df in dfs:
		# binned_dfs.append(bin_dataframe(dfm 
	
	print "plotting cumulative sum"
	fig = plot_cum_sum(dfs, args.groupnames)
	save_figure(fig, args.out+'.cumsum', args.figformat)
	
	print "plotting radial and 2D intensity"
	fig = plot_raw_data(dfs, args.groupnames, args.out+'.hist2d', rbins=1000, ibins=1000)
	save_figure(fig, args.out, args.figformat)
#end compare_raw_data()

def coloc_raw_data(args):
	dfs = fetch_raw_data(args)
	
	print "plotting colocalization results"
	fig = plot_colocalization_results(dfs, args.out, dfs[0].columns[1], dfs[0].columns[2], t1=0, t2=0.4, labels=args.groupnames, ibins=1000)
	save_figure(fig, args.out+'.coloc', args.figformat)
	
	print "plotting cumulative sum"
	fig = plot_cum_sum(dfs, args.groupnames)
	save_figure(fig, args.out+'.cumsum', args.figformat)
	
	print "plotting radial and 2D intensity"
	fig = plot_raw_data(dfs, rbins=1000, ibins=1000)
	save_figure(fig, args.out, args.figformat)
#end coloc_raw_data()

def plot_cum_sum(dfs, labels):
	fig = plt.figure(1, figsize=(12, (4 * len(dfs))), dpi=600)
	plt.subplots_adjust(hspace=0.6, wspace=0.5)
	gridsize = (1,1)
	
	ax = plt.subplot2grid(gridsize, (0,0))
	ip.plot_radial_cumulative_sum(ax, dfs, dfs[0].columns[0], dfs[0].columns[1], labels=labels)
	return fig
#end plot_cum_sum()

def plot_raw_data(dfs, labels, hdatabasefn, rbins=500, ibins=500):
	fig = plt.figure(1, figsize=(12, (4 * len(dfs))), dpi=600)
	plt.subplots_adjust(hspace=0.6, wspace=0.5)
	df_cols = len(dfs[0].columns)
	combinations = list(itertools.combinations([dfs[0].columns[i] for i in range(1, len(dfs[0].columns))], 2))
	gridsize = (len(dfs), (len(dfs[0].columns)-1) + len(combinations))
	
	for d in range(len(dfs)):
		df = dfs[d]
		for i in range(1, df_cols):
			ax = plt.subplot2grid(gridsize, (d,i-1))
			hist_data = ip.plot_radial_intensity(ax, df, df.columns[0], df.columns[i], radial_bins=rbins, intensity_bins=ibins)
			write_2D_hist_data(hist_data['hist2d'], '{}.{}.tsv'.format(hdatabasefn, labels[d]), labels[d], 1)
			
		for j in range(len(combinations)):
			ax = plt.subplot2grid(gridsize, (d,df_cols+j-1))
			ip.plot_2D_intensity(ax, df, combinations[j][0], combinations[j][1], bins=ibins)
		
	return fig
#end plot_raw_data()

def write_2D_hist_data(hist_data, filename, genotype, expid):
	hist = hist_data[0]
	xbins = hist_data[1]
	ybins = hist_data[2]
	
	template = '{geno}\t{expid}\t{rad}\t{inten}\t{cnt}\n'
	
	with open(filename, 'w+') as of:
		of.write('Geno\tExp\tRad\tInten\tCount\n')
		for i in range(len(ybins)-1):
			for j in range(len(xbins)-1):
				of.write(template.format(geno=genotype, expid=expid, rad=xbins[j], inten=ybins[i], cnt=hist[i,j]))

#end write_2D_hist_data()

def plot_colocalization_results(dfs, dest_base, c1, c2, t1, t2, labels, ibins=500):
	fig = plt.figure(1, figsize=((4 * len(dfs)), 12), dpi=600)
	plt.subplots_adjust(hspace=0.6, wspace=0.5)
	gridsize = (3, len(dfs))
	half_width = (len(dfs) / 2)
	
	thold_dfs = [df[(df[c1]>=t1) & (df[c2]>=t2)] for df in dfs]
	for d in range(len(dfs)):
		df = thold_dfs[d]
		ax_hist2D = plt.subplot2grid(gridsize, (0,d))
		ip.plot_2D_intensity(ax_hist2D, df, c1, c2, bins=ibins)
		
	ax_bar = plt.subplot2grid(gridsize, (1,0), rowspan=2, colspan=half_width)
	ax_bar.set_title('Mean {} Intensity in voxels with {} >= {}'.format(c1, c2, t1))
	ip.plot_bar_intensities(ax_bar, thold_dfs, c1, labels)
	
	ax_vio = plt.subplot2grid(gridsize, (1,half_width), rowspan=2, colspan=half_width)
	ax_vio.set_title('{} Intensity in voxels with {} >= {}'.format(c1, c2, t1))
	ip.plot_violin_intensities(ax_vio, thold_dfs, c1)
	
	return fig
#end plot_colocalization_results



def compare_binned_data(args):
	dfs = []
	for i in range(len(args.files)):
		sys.stderr.write('Reading file {}....\n'.format(args.files[i]))
		dfs.append({
			'i': i,
			'n': get_label_from_path(args.files[i]),
			'c': ip.get_color(i),
			'd': pd.read_csv(args.files[i], sep='\t', comment='#', skip_blank_lines=True, index_col='Radial')
		})
	
	fig = plt.figure(1, figsize=(11, 8), dpi=600)
	plt.subplots_adjust(hspace=0.6, wspace=0.5)
	gridsize = (3, 20)
	x = dfs[0]['d'].index.values
	
	
	#begin scatter/line plot
	ax_scat = plt.subplot2grid(gridsize, (0,0), colspan=gridsize[1], rowspan=2)
	i=0
	for df in dfs:
		y = ttstats.summarize_data(df['d'], method=args.stat, axis=1)
		computed_error = ttstats.compute_error(df['d'], method='sem', axis=1)
		ax_scat.plot(x, y, label=df['n'], color=df['c'])
		ax_scat.fill_between(x, y, (y + computed_error[0]), facecolor=df['c'], edgecolor='none', alpha=0.2)
		ax_scat.fill_between(x, y, (y - computed_error[1]), facecolor=df['c'], edgecolor='none', alpha=0.2)
		i += 1
	ax_scat.set_xlim(0, 1)
	legend = ax_scat.legend(loc='upper left')
	#end scatter/line plot
	
	
	#begin violin plot
	for i in range(len(x)):
		ax_vio = plt.subplot2grid(gridsize, (2,i), colspan=1, rowspan=2)
		summary = [df['d'].iloc[i].values for df in dfs]
		#print summary
		vparts = ax_vio.violinplot(summary, positions=[df['i'] for df in dfs], showmeans=False, showmedians=False, showextrema=False)
		
		#change the filled area of the violin
		for df in dfs:
			vparts['bodies'][df['i']].set_facecolor(df['c'])
			vparts['bodies'][df['i']].set_alpha(1)
			
			quartile1, medians, quartile3 = np.percentile(summary[df['i']], [25, 50, 75])
			whiskers = ip._adjacent_values(sorted(summary[df['i']]), quartile1, quartile3)
			whiskersMin, whiskersMax = whiskers[0], whiskers[1]
			ax_vio.scatter([df['i']], medians, marker='o', color='white', s=0.1, zorder=3)
			ax_vio.vlines([df['i']], quartile1, quartile3, color='k', linestyle='-', lw=1, alpha=0.5)
			ax_vio.vlines([df['i']], whiskersMin, whiskersMax, color='k', linestyle='-', lw=0.25, alpha=0.5)
			
		ax_vio.set_title("{:0.2f}".format(x[i]), rotation=90, verticalalignment='bottom')
		ax_vio.set_ylim(0,0.8)
		ax_vio.set_xticklabels([])
		ax_vio.set_xticks([])
		if i > 0:
			ax_vio.set_yticklabels([])
			ax_vio.set_yticks([])
	#end violin plot
	
	#finally save the figure
	save_figure(fig, args.out, args.figformat)
#end compare_data()






#####################
# UTILITY FUNCTIONS #
#####################
def save_figure(fig, basename, fig_format):
	savename = "{}.{}".format(basename, fig_format)
	sys.stderr.write('Saving Figure.....\n')
	fig.savefig(savename, dpi=600)
	plt.close(fig)
	sys.stderr.write(' => See %s for results.\n' % (savename,))
#end save_figure()

def get_label_from_path(path):
	return os.path.splitext(os.path.basename(path))[0]
#end et_label_from_path()

def normalize_data(data, columns, ignore_outliers=True):
	for c in columns:
		if ignore_outliers:
			datamax = ttstats.max_exclude_outliers(data[c])
		else:
			datamax = data[c].max()
		data[c] /= datamax
	return data
#end normalize_data()













if __name__ == "__main__":
	main()








