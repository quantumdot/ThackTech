import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
from ThackTech.Plotting import stats as ttstats
from matplotlib.colors import LogNorm
from xlrd.xlsx import cell_name_to_rowx_colx



color_cycle = ['navy','blue','dodgerblue','firebrick','tomato','lightsalmon','k']
def get_color(i):
    '''Gets the next color in the color cycle based on index i
    
    Parameters:
        i: int index of the color to get
        
    Returns:
        a valid color for use in matplotlib
    '''
    return color_cycle[i % len(color_cycle)]
#end get_color()


def plot_radial_intensity(ax, df, col_rad, col_int, radial_bins=1000, intensity_bins=1000, **kwargs):
    ''' Wrapper method for plot_2D_intensity() that sets a bunch of 
        useful defaults for plotting radial position (x) vs intensity (y).
        
        radial_bins and intensity_bins are proxies for xbins and ybins, respectively.
        Other arguments to plot_2D_intensity() can be overridden via **kwargs
    '''
    args = {
        'xbins': radial_bins, 
        'ybins': intensity_bins,
        'xrange': (0,1),
        'yrange': (0,1),
        'title': col_int,
        'xlabel': 'Radial Position',
        'ylabel': 'Relative Intensity',
        'linreg': True,
        'avg': True
    }
    args.update(kwargs)
    return plot_2D_intensity(ax, df, col_rad, col_int, **args)
#end plot_radial_intensity()
    
def plot_2D_intensity(ax, df, xcol, ycol, xbins=1000, ybins=1000, xrange=(0,1), yrange=(0,1), title=None, xlabel=None, ylabel=None, linreg=True, corr=True, avg=True, **kwargs):
    '''Plots a 2D histogram on the supplied axis
    
    Parameters:
        ax: axis to plot 2D histogram on
        df: dataframe containing data to plot
        xcol: column in dataframe that represents x-axis data
        ycol: column in dataframe that represents y-axis data
        xbins: number of bins to use on the x-axis data
        ybins: number of bins to use on the y-axis data
        xrange: left and right-most edges of bins in the x-axis
        yrange: left and right-most edges of bins in the y-axis
        title: title for the plot. If None, constructed from x and y column labels
        xlabel: label for the x-axis of the plot. If None, constructed from the x column label
        ylabel: label for the y-axis of the plot. If None, constructed from the y column label
        linreg: if True, plot a linear regression of the data
        avg: if True, plot a binned mean across the x-axis
        **kwargs: additional arguments to pass to hist2D method
        
    Returns:
        Dict of the following items:
            hist2d: Tuple contianing the following items:
                - counts: 2D array of histogram counts
                - xedges: 1D array of x edges
                - yedges: 1D array of y edges
                - Image: 2D array of image data
            linreg: Present if linreg is True. Tuple containing the following items:
                - m: slope of the linear regression
                - b: y-intercept of the linear regression
            corr: Present if corr is True. Dictionary containing the following items:
                - spearman: Tuple containing spearman correlation coefficient and pvalue
                - pearson: Tuple containing pearson correlation coefficient and pvalue
            avg: Present if avg is True. Tuple containing the following items:
                - statistic: values of the mean statistic for each bin
                - bin_edges: bin edges
                - binnumber: an integer for each observation that represents the bin in which this observation falls. Array has the same length as values.
        
    '''
    return_values = {}
    #plot the 2D-histogram with colorbar
    ax_H = ax.hist2d(df[xcol], df[ycol], bins=[xbins, ybins], norm=LogNorm(), range=[xrange,yrange], **kwargs)#, alpha=0.1, s=10, linewidths=0)
    plt.colorbar(ax_H[3], ax=ax, label="Frequency")
    return_values['hist2d'] = ax_H
    
    if linreg:
        #make the linear regression
        m, b = np.polyfit(df[xcol], df[ycol], 1)
        ax.plot(df[xcol], m * df[xcol] + b, '-', c='k')
        return_values['linreg'] = (m, b)
        
    if corr:
        return_values['corr'] = {
			"spearman": stats.spearmanr(df[xcol], df[ycol]),
			"pearson": stats.pearsonr(df[xcol], df[ycol])
		}
     
    if avg:
        #plot the mean of y along the x
        data_stats, bin_edges, binnumber = stats.binned_statistic(df[xcol], df[ycol], statistic='mean', bins=xbins, range=xrange)
        ax.plot(bin_edges[1:], data_stats, ':', color='k')
        return_values['avg'] = (data_stats, bin_edges, binnumber)
    
    #make it pretty and informative
    ax.set_title(xcol + ' vs. ' + ycol) if (title is None) else ax.set_title(title)
    ax.set_xlabel(xcol + ' Intensity') if xlabel is None else ax.set_xlabel(xlabel)
    ax.set_ylabel(xcol + ' Intensity') if ylabel is None else ax.set_ylabel(ylabel)
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.set_aspect('equal')
    
    return return_values
#end plot_2D_intensity()

def plot_violin_intensities(ax, dfs, col):
    summary = [df[col].values for df in dfs]
    vparts = ax.violinplot(summary, positions=range(len(dfs)), showmeans=False, showmedians=False, showextrema=False)
    
    #change the filled area of the violin
    for i in range(len(dfs)):
        vparts['bodies'][i].set_facecolor(color_cycle[i % len(color_cycle)])
        vparts['bodies'][i].set_alpha(1)
        
        quartile1, medians, quartile3 = np.percentile(summary[i], [25, 50, 75])
        whiskers = _adjacent_values(sorted(summary[i]), quartile1, quartile3)
        whiskersMin, whiskersMax = whiskers[0], whiskers[1]
        ax.scatter([i], medians, marker='o', color='white', s=0.1, zorder=3)
        ax.vlines([i], quartile1, quartile3, color='k', linestyle='-', lw=1, alpha=0.5)
        ax.vlines([i], whiskersMin, whiskersMax, color='k', linestyle='-', lw=0.25, alpha=0.5)
    
    #make it pretty and informative
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylim(0,1.0)
    #if i > 0:
    #    ax.set_yticklabels([])
    #    ax.set_yticks([])
#end plot_violin_intensities()

def _adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value
#end adjacent_values()

def plot_bar_intensities(ax, dfs, col, labels):
    for i in range(len(dfs)):
        df = dfs[i]
        summary = df[col].mean()
        error = df[col].sem()
        ax.bar([i], [summary], yerr=[error], label=labels[i], color=color_cycle[i % len(color_cycle)])
    ax.legend()
    ax.set_xticklabels([])
    ax.set_xticks([])
#end plot_bar_intensities

def plot_line_intensities(ax, dfs, colx, coly, labels, summary_method, error_method, axis=1):
    for i in range(len(dfs)):
        df = dfs[i]
        color = color_cycle[i % len(color_cycle)]
        y = ttstats.summarize_data(df[coly], method=summary_method, axis=axis)
        x = df[colx]
        computed_error = ttstats.compute_error(df[coly], method=error_method, axis=axis)
        ax.plot(x, y, label=labels[i], color=color)
        ax.fill_between(x, y, (y + computed_error[0]), facecolor=color, edgecolor='none', alpha=0.2)
        ax.fill_between(x, y, (y - computed_error[1]), facecolor=color, edgecolor='none', alpha=0.2)

    #ax_scat.set_xlim(0, 1)
    legend = ax.legend(loc='best')
#end plot_line_intensities()

def plot_radial_cumulative_sum(ax, dfs, colx, coly, labels):
    
    for i in range(len(dfs)):
        print "plotting cumulative radial sum for {}".format(labels[i]) 
        df = dfs[i]
        color = color_cycle[i % len(color_cycle)]
        df_s = df.sort_values(colx)
        df_s = df_s.where(df_s[colx] <= 1.0)
        
        df_s['cum_sum'] = df_s[coly].cumsum()
        df_s['cum_perc'] = df_s['cum_sum']/df_s[coly].sum()
        #print df_s
        #print df_s[coly].sum()
        
        print "radial\t%sig"
        for j in np.arange(0.0, 1.0, 0.05, dtype=float):
            print "{}\t{}".format(j, df_s.where(df_s[colx]>=j)['cum_perc'].min())
            
        print "%sig\tradial"
        for j in np.arange(0.0, 1.0, 0.05, dtype=float):
            print "{}\t{}".format(j, df_s.where(df_s['cum_perc']>=j)[colx].min())
        
        ax.plot(df_s[colx], df_s['cum_perc'], label=labels[i], color=color)

    ax.set_xlabel('Radial Position')
    ax.set_ylabel('Cumulative Percent Sum Intensity')
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    legend = ax.legend(loc='best')
#end plot_line_intensities()


