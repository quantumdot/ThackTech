#Profiler Usage

## Parameters

### Input Data Options:

| Option | Description | Default |
|-------------------|------------------------------------------------------------------------------------------------------------------------------|-------|
| `-b BED, --bed BED` | Bed files to profile against. | None |
| `-s SIG, --sig SIG` | Signal files to profile. BigWig, BAM, BED, BigBed files are supported. All but bigwig will be treated by counting intervals. | None |
| `-i INP, --inp INP` | Input files to be used for normalization. | None |
| `--ignoreinput` | Do not use input data. | False |

### Labeling Options:

| Option | Description | Default |
|-------------------------------|------------------------------------------------------------------------------------------------------------------------------|---------|
| `-sl SLABEL, --slabel SLABEL` | Signal labels for each plot. Specify in the same order as the signal files. If not supplied, file basename will be used. | None |
| `-il ILABEL, --ilabel ILABEL` | Interval labels for each plot. Specify in the same order as the interval files. If not supplied, file basename will be used. | None |
| `--title TITLE` | Title for entire plot | "" |
| `--bedscorelabel BEDSCORELABEL` | Label for plotting the score column of intervals. | "Score" |
| `--xnumticks XNUMTICKS` | Number of tick marks to use in the x-axis. | 3 |
| `--ynumticks YNUMTICKS` | Number of tick marks to use in the y-axis. | 4 |
| `--xlabelrot XLABELROT` | Angle, in degrees, for x axis label rotation. | 0 |
| `--ylabelrot YLABELROT` | Angle, in degrees, for y axis label rotation. | 90 |


### Profiling Options:

| Option | Description | Default |
|------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|
| `--nochromfilter` | Do not filter bed file for common chromosomes. By default profiling only occurs on chromosomes common to all files (useful for ignoring random or Un* chromosomes). | False |
| `--chrignore CHRIGNORE` | Regex for chromosome strings to filter out of intervals. | None |
| `--up UP` | Span upstream (in bp) from selected intervals alignement to profile. | 1000 |
| `--down DOWN` | Span downstream (in bp) from selected intervals alignement to profile. | 1000 |
| `--res RES` | Profiling resolution in bp. | 10 |
| `--dir` | If set, the direction (+/-) [strand orientation] is considered in profiling. Strand information should be present in bed file used, and if not present is assumed to be +. | False |
| `--align {left,center,right,scale}` | Method used to align intervals. Left aligns intervals along the 5' end, while right aligns intervals aling the 3' end (assuming intervals provide strand information and --dir is specified, otherwise + strand is assumed). Center aligns intervals along the center center point of intervals. Scale will scale all intervals to the same apparent width specified by --scaleregionsize. | center |
| `--scaleregionsize SCALEREGIONSIZE` | Size of the scaled region (i.e. gene-body) in bp. This option is only useful for "--align scale" option. | 3000 |
| `--nan {zero,ignore}` | How to handle missing or NaN values. | ignore |
| `--collectionmethod {get_as_array,ucsc_summarize,summarize}` | Method for collecting signal data. | get_as_array |
  

### Scaling Options:

| Option | Description | Default |
|-----------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `--scalegroups SCALEGROUPS` | Groups of plots to share color/y-axis scales. If not specified, all plots will be constructed with the same scale. Parameter should be specified as comma- separated lists of 0-based offsets of samples, and groups separated with a semicolon. Ex: 0;1,2;3,4 results in sample 0 plotted with independent scale, 1 and 2 sharing scale, and 3 and 4 sharing scale. If specified, but parameter omits samples, then the omitted samples will each be scaled independently. | None |
| `--scalebedgroups SCALEBEDGROUPS` | Groups of plots to share color/y-axis scales. If not specified, all plots will be constructed with the same scale. Parameter should be specified as comma- separated lists of 0-based offsets of samples, and groups separated with a semicolon. Ex: 0;1,2;3,4 results in sample 0 plotted with independent scale, 1 and 2 sharing scale, and 3 and 4 sharing scale. If specified, but parameter omits samples, then the omitted samples will each be scaled independently. | None |
| `--saturatemin SATURATEMIN` | In the heatmap plot, saturate the <--saturatemin> percent bottom values. | 0.01 |
| `--saturatemax SATURATEMAX` | In the heatmap plot, saturate the <--saturatemax> percent top values. | 0.01 |
| `--coefficients [COEFFICIENTS [COEFFICIENTS ...]]` | Coefficients to multiply signals by, one value for each signal submitted. | None |
| `--normalizationmethod {ratio,log2,reciprocal_ratio,subtract,add,mean}` | Method to use for normalizing signal by input. | log2 |

### Clustering/Sorting Options:

| Option | Description | Default |
|----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `--sort SORT` | sample index (0-based) to use for sorting. If not specified than the order of the bed file is used. Mutually exclusive with --kmeans. | None |
| `--sortmethod {mean,median,max,min,sum}` | Method used for sorting. | "mean" |
| `--sortrange SORTRANGE` | Range of the profiles (in relative bp) to be used in the sorting operation. Specify in the format "start:stop". Default is to use the entire range. | None |
| `--kmeans` | If set, Perform K-means clustering on the data. Mutually exclusive with --sort. | False |
| `--k K` | Number of clusters, k, to fit data to when performing K-means clustering. | None |
| `--autok` | Optimize the number of clusters, k, in the dataset. Mutually exclusixe with --k. | False |
| `--ksamples KSAMPLES` | Comma-separated list of 0-based offsets of samples to use for K-means clustering. Use 'all' to cluster on all samples. | "all" |
| `--summarymethod` | Method used for producing summary (avg) plot. | "mean" |
| `--hline` | Draw horizontal lines to delineate clusters. | False |
| `--hlineweight HLINEWEIGHT` | Horizontal line weight. | 0.5 |
| `--hlinestyle HLINESTYLE` | Line style for the horizontal lines. See matplotlib axhline documentation for valid values. | "-" |

### Output Options:

| Option | Description | Default |
|-------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `--name NAME` | Base name for the plot output. | "" |
| `--format {pdf,png,svg}` | Format to output the final figure. | "pdf" |
| `--plot {avg,violin,heat,kavg,bedscores,truebedscores}` | Types of plots to produce. Supply multiple options to produce hybrid plots. avg will produce average profiles, heat will produce heatmaps, kavg will produce average profiles for each class determined by k-means clustering (only available when used with --plot heat AND --kmeans). bedscores will plot the score column for intervals simply, while truebedscores will plot the score column for intervals more accuratly (showing true genomic context of element; mutually exclusive with --plot bedscores). | None |
| `--dpi DPI` | DPI resolution of the saved figure. | 600 |
| `--width WIDTH` | Width (in inches) of the figure. | 8.0 |
| `--height HEIGHT` | Height (in inches) of the figure. | 6.0 |
| `--avgplotrows AVGPLOTROWS` | Number of rows to use for average profile plots. | 1 |
| `--heatplotrows HEATPLOTROWS` | Number of rows to use for heatmap plots. | 3 |
| `--rotate` | By default plots will be arranged with signals going across and intervals going down. --rotate will change the orientation so that intervals go across and signals going down. | False |
| `--dumpsummary` | Write summary data to a text file. See --summaryrange and --summarymethod for details on affecting the output generated by this switch. | False |
| `--summaryrange SUMMARYRANGE` | Range(s) of the profiles (in relative bp) to be used in the sorting operation. Specify in the format "start:stop". Default is to use the entire range. Write multiple ranges by separating ranges with a semicolon. | None |

### General Plotting Options:

| Option | Description | Default |
|---------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| `--vline` | Draw a vertical line at the "zero" point. | False |
| `--vlineweight VLINEWEIGHT` | Vertical line weight. | 0.5 |
| `--vlinestyle VLINESTYLE` | Line style for the vertical line. See matplotlib axvline documentation for valid values. | ":" |
| `--fontsize FONTSIZE` | Font size used in plots. | 8 |
| `--legendfontsize LEGENDFONTSIZE` | Font size used in plot legend. | 6 |
| `--showci` | Plot confidence interval for average plot. | False |
| `--ciwidth CIWIDTH` | Confidence interval width for average plot. sem for Standard Error of the mean, std for Standard deviation, or float for percent CI (i.e. 0.95). | sem |
| `--genome GENOME` | UCSC reference genome that intervals and signal are computed with (i.e. mm9, hg19, etc.). This is really only required when using --plot truebedscores. | None |
| `--colors COLORS [COLORS ...]` | Colors to cycle through for individual samples on aggregate plots. | ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'firebrick', 'darkolivegreen', 'navy', 'palevioletred'] |
| `--linewidth LINEWIDTH` | Width of line used for aggregate plots (i.e. average profile). | 0.5 |

### Resource Options:

| Option | Description | Default |
|----------------------------|--------------------------------------|-------------------|
| `--cache` | If set, profiles are dumped to disk. | False |
| `--cachedir CACHEDIR` | Location to place cache files. | ".profiler_cache" |
| `-p CPUS, --processors CPUS` | Number of processors to use. | 4 |


## Examples:
coming soon....
