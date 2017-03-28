#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import os
import sys
import argparse
import pybedtools
import matplotlib.pyplot as plt
import numpy
from scipy import stats
from matplotlib.ticker import NullFormatter


class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
        
    def __del__(self):
        self.release()
        
    def release(self):
        sys.stdout = self.stdout
        self.file.close()
        
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

class AutoAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not values:
            setattr(namespace, self.dest, ['auto'])
        else:
            setattr(namespace, self.dest, values)


class DescriptiveStats(object):
    def __init__(self, data):
        self.data = data
        self.compute()

    def compute(self):
        self.count  = self.data.size
        self.min    = self.data.min()
        self.max    = self.data.max()
        self.mean   = self.data.mean()
        self.median = numpy.median(self.data)
        self.mode   = stats.mode(self.data)[0][0]
        self.stdev  = self.data.std()
    
    def render(self):
        buffer  = ""
        buffer += "Count     = "+str(self.count)+"\n"
        buffer += "Min       = "+str(self.min)+"\n"
        buffer += "Max       = "+str(self.max)+"\n"
        buffer += "Mean      = "+str(self.mean)+"\n"
        buffer += "Median    = "+str(self.median)+"\n"
        buffer += "Mode      = "+str(self.mode)+"\n"
        buffer += "Std. Dev. = "+str(self.stdev)+"\n"
        return buffer        



class MetricContainer(object):
    def __init__(self, name, data):
        self.name = name
        self.data = data
        self.stats = DescriptiveStats(self.data)
    
    def render_stats(self):
        buffer  = "Statistics for "+self.name+" Distribution:\n"
        buffer += "----------------------------------\n"
        buffer += self.stats.render()
        buffer += "----------------------------------\n\n"
        return buffer

    def plot(self, ax, bins, log, limit):
        ax.hist(self.data, bins=bins, color='b')
        plt.title("Interval "+self.name+" Distribution")
        plt.xlabel(self.name)
        if log:
            ax.set_yscale('log', nonposy='clip')
            plt.ylabel("log(Frequency)")
        else:
            plt.ylabel("Frequency")

        if limit is not None:
            ax.set_xlim(right=limit)


class cxt.sampleContainer(object):
    def __init__(self, name):
        self.name = name
        self.metrics = {}
    
    def addMetric(self, name, data):
        self.metrics[name] = MetricContainer(name, data)



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', nargs='+', help="BED file(s) to calculate statistics.")
    parser.add_argument('--bins', default=1000, type=int, help="Number of bins in the histogram. (default = 1000)")
    parser.add_argument('--format', default='png', help="Output format, one of [png, pdf, ps, eps, svg]. (defualt = png)")
    parser.add_argument('--xkcd', action='store_true', help="Create the plots in the xkcd style")
    parser.add_argument('--D2', action='store_true', help="Produce a 2D plot using bed scores and lengths.")
    parser.add_argument('--log', action='store_true', help="Use log scale in the frequency axis of histograms.")
    parser.add_argument('--squeeze', action=AutoAction, default=[None], type=str, nargs='*', help="Force the axis max limits to be the same amongst all cxt.samples. If not specified then each cxt.sample will have independent ranges. If specified with no argument or 'auto', the range will be deduced from all the max of all cxt.samples for each distribution. You may also specify a number that will be used for the range. Multiple values affect each distribution in [alphebitical] order. You may mix 'auto' and numbers.")
    args = parser.parse_args()

    #print args
    #exit()

    cxt.samples = []
    for b in args.bed:
        cxt.samples.append(parse_bed(b))

    metrics = sorted(cxt.samples[0].metrics.keys())
    limits = {}
    mi = 0
    for mk in metrics:
        if (args.squeeze[mi] is None) or (args.squeeze[mi] == 'none'):
            #no squeezing is done
            limits[mk] = None
        elif args.squeeze[mi] == 'auto':
            #automatic squeezing
            limits[mk] = find_auto_limit(cxt.samples, mk)
        else:
            #we should have gotten a number from the user
            limits[mk] = float(args.squeeze[mi])
        mi = clamp(mi+1, 0, len(args.squeeze)-1)
        
    for s in cxt.samples:
        process_cxt.sample(s, args, limits)
#end main()

def clamp(val, minval, maxval):
    if val < minval: return minval
    if val > maxval: return maxval
    return val
#end clamp()

def find_auto_limit(cxt.samples, metric):
    limit = 0
    for s in cxt.samples:
        limit = max(limit, s.metrics[metric].stats.max)
    return limit
#end find_auto_limit()

def parse_bed(bedfile):
    sys.stdout.write("Reading "+bedfile+".....\n")
    bed = pybedtools.BedTool(bedfile)
    savename = os.path.splitext(os.path.basename(bedfile))[0]
    cxt.sample = cxt.sampleContainer(savename)

    #initialize some vars...
    lengths = numpy.empty(bed.count(), dtype=int)
    scores = numpy.empty(bed.count(), dtype=float)
    i = 0
    for interval in bed:
        lengths[i] = interval.length
        scores[i] = interval.score
        i += 1
    cxt.sample.addMetric('Length', lengths)
    cxt.sample.addMetric('Score', scores)
    return cxt.sample
#end parse_bed()


def process_cxt.sample(cxt.sample, args, limits):

    sys.stdout = Tee(cxt.sample.name+'.stats.txt', 'w')
    sys.stdout.write("\nProcessing "+cxt.sample.name+".bed....\n\n")

    #allow xkcd style graphs
    if args.xkcd:
        plt.xkcd()

    plt.figure(1, figsize=(8, 11))
    plt.subplots_adjust(hspace=1)

    gridsize = (1,1)
    if args.D2:
        gridsize = (4,1)

    lengths = cxt.sample.metrics['Length']
    scores  = cxt.sample.metrics['Score']

    #plot the lengths    
    ax1 = plt.subplot2grid(gridsize, (0,0))
    lengths.plot(ax1, args.bins, args.log, limits[lengths.name])
    sys.stdout.write(lengths.render_stats())

    if args.D2:
        ax2 = plt.subplot2grid(gridsize, (1,0))
        scores.plot(ax2, args.bins, args.log, limits[scores.name])
        sys.stdout.write(scores.render_stats())

            
        #2D scatter plot
        ax3 = plt.subplot2grid(gridsize, (2,0), rowspan=2)
        ax3.scatter(lengths.data, scores.data, alpha=0.3, s=10, linewidths=0)
        ax3.set_xlim(left=0)
        ax3.set_ylim(bottom=0)
        if args.squeeze:
            ax3.set_xlim(right=limits[lengths.name])
            ax3.set_ylim(top=limits[scores.name])
        
        plt.title("Interval Length vs. Interval Scores")
        plt.xlabel("Interval Length")

    imgname = cxt.sample.name+'.stats'+'.'+args.format
    plt.savefig(imgname)
    sys.stdout.write("See graphical output at: "+imgname+"\n\n")
    sys.stdout.release() #release the logger!
#end process_bed()




if __name__ == "__main__":
    main()
