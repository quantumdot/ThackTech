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


def plot_radial_intensity(ax, df, col_rad, col_int, radial_bins=1000, intensity_bins=1000):
    #plot the 2D-histogram with colorbar
    ax_H = ax.hist2d(df[col_rad], df[col_int], bins=[radial_bins, intensity_bins], norm=LogNorm())
    plt.colorbar(ax_H[3], ax=ax, label="Frequency")
    #make the linear regression
    ax_m, ax_b = np.polyfit(df[col_rad], df[col_int], 1)
    ax.plot(df[col_rad], ax_m * df[col_rad] + ax_b, '-', c='k')
    #plot the mean of y along the x
    data_stats, bin_edges, binnumber = stats.binned_statistic(df[col_rad], df[col_int], statistic='mean', bins=radial_bins, range=[(0, 1)])
    ax.plot(bin_edges[1:], data_stats, ':', color='k')
    #make it pretty and informative
    ax.set_title(col_int)
    ax.set_xlabel('Radial Position')
    ax.set_ylabel('Relative Intensity')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
#end plot_radial_intensity()
    
def plot_2D_intensity(ax, df, colx, coly, bins=1000):
    #plot the 2D-histogram with colorbar
    ax_H = ax.hist2d(df[colx], df[coly], bins=bins, norm=LogNorm(), range=[[0,1],[0,1]])#, alpha=0.1, s=10, linewidths=0)
    plt.colorbar(ax_H[3], ax=ax, label="Frequency")
    #make the linear regression
    m, b = np.polyfit(df[colx], df[coly], 1)
    ax.plot(df[colx], m * df[colx] + b, '-', c='k')
    #make it pretty and informative
    ax.set_title(colx + ' vs. ' + coly)
    ax.set_xlabel(colx + ' Intensity')
    ax.set_ylabel(coly + ' Intensity')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
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
        df = dfs[i]
        color = color_cycle[i % len(color_cycle)]
        df_s = df.sort_values(colx)
        df_s['norm_x'] = df_s[coly] / df_s[coly].max()
        df_s['cum_sum'] = df_s[coly].cumsum()
        df_s['cum_perc'] = df_s['cum_sum']/df_s[coly].sum()
        print df_s
        print df_s[coly].sum()
        
        
        ax.plot(df_s['norm_x'], df_s['cum_perc'], label=labels[i], color=color)

    #ax.set_xlim(0, 1.0)
    #ax.set_ylim(0, 1.0)
    legend = ax.legend(loc='best')
#end plot_line_intensities()



