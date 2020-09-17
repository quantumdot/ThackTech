#!/usr/bin/python
"""
Given an arbitrary number of files, computes the combinatorial intersections of these files.

The values in the diagram assume:

    * unstranded intersections
    * no features that are nested inside larger features

    *this isnt absolutly perfect. see this discussion on why the numbers wont always add up!
    --https://github.com/daler/pybedtools/issues/45

"""

import argparse
import sys
import os
import pybedtools
import itertools
import csv
import json
import subprocess
import math
from ThackTech.Processes import ProgressBar


def main():

    op = argparse.ArgumentParser(prog=sys.argv[0])
    op.add_argument('bed', nargs='+',  help='Bed files to compute intersections')
    op.add_argument('-c', action='append', help='Optional colors for circles. Specify as Hex color or in "rgb(r,g,b)" format.')
    op.add_argument('-l', action='append', help='Labels for bed files. Specify in the same order as the BEDs. If not supplied, BED basename will be used.')
    #op.add_argument('--size', default='300x300', help='Optional size of PNG, in pixels.  Default is "%(default)s"')
    op.add_argument('-o', default='out', help='Output file basename to save as')
    op.add_argument('-m', default='classic', choices=['classic', 'edwards'], help="Sets the display mode for the chart output.")
    op.add_argument('--pdf', action='store_true', help="Convert the HTML file automagically to PDF. Requires the wkhtmltopdf package.")
    options = op.parse_args()

    
    beds = [None]*len(options.bed)
    labels = [None]*len(options.bed)
    bi = 0
    for b in options.bed:
        labels[bi] = options.l[bi] if (options.l is not None and bi < len(options.l)) else os.path.basename(b)
        beds[bi] = pybedtools.BedTool(b)
        bi += 1
    
    progress = ProgressBar(math.pow(len(beds), 2) - 1, "Computing Intersections....", barlength=50)
    progresscounter = 0

    #print labels    
    #compute the matrix
    letters = ['A', 'B', 'C', 'D', 'E', 'F']
    data = {}
    chartdata = { 'name': {}, 'values': {} }
    for i in range(1, len(beds)+1):
        combos = itertools.combinations(range(len(beds)), i)
        for tup in combos:
            progress.update(progresscounter, str(tup))
            #print "Processing Tuple: "+str(tup)
            l = labels[tup[0]]
            cl = letters[tup[0]]
            d = beds[tup[0]]
            eli = 0
            #record the ANDS 
            for el in tup:
                if eli > 0:
                    l += " + "+labels[el]
                    cl += letters[el]
                    d += beds[el]
                eli += 1
            #compute the NOTS
            compliment = getComplement(range(len(beds)), tup)
            for c in compliment:
                d -= beds[c]
            data[l] = d.count()
            chartdata['values'][cl] = data[l]
            
            #report Progress
            progresscounter += 1
            progress.update(progresscounter)
    
    for i in range(len(beds)):
        chartdata['name'][letters[i]] = labels[i]

    chartsetup = {}
    chartsetup['series'] = [chartdata]
    chartsetup['displayMode'] = options.m
    chartsetup['displayStat'] = True
    if options.c:
        chartsetup['colors'] = options.c
    
    chartbasename = options.o+'.'+options.m
    writeChart(chartsetup, chartbasename)

    if options.pdf:
        print "Generating PDF....."
        subprocess.call("wkhtmltopdf --quiet "+chartbasename+".html "+chartbasename+".pdf", shell=True) 

    writer = csv.writer(open(options.o+'.csv', 'wb'))
    for key, value in data.items():
        writer.writerow([key, value])


def getComplement(items, selected):
    complement = []
    for item in items:
        if item not in selected:
            complement.append(item)
    return complement

    

def writeChart(json_data, filename):
    js_src = [
        'http://jvenn.toulouse.inra.fr/app/js/jquery.min.js',
        'http://jvenn.toulouse.inra.fr/app/js/bootstrap-colorpicker.min.js',
        'http://jvenn.toulouse.inra.fr/app/js/canvas2svg.js',
        'http://jvenn.toulouse.inra.fr/app/js/jvenn.min.js'
    ]
    css = [
        'http://jvenn.toulouse.inra.fr/app/css/bootstrap.css',
        'body { padding-top: 60px; padding-bottom: 40px; }',
        'http://jvenn.toulouse.inra.fr/app/css/prettify.css',
        'http://jvenn.toulouse.inra.fr/app/css/bootstrap-responsive.css'
    ]
    cpyright = '<p>Venn Diagram Rendering: Copyright &copy; 2014, INRA | Designed by <a href="http://bioinfo.genotoul.fr" target="_blank">GenoToul Bioinfo</a> and<a href="http://sigenae.org" target="_blank">Sigenae</a> teams.</p>'
    cpyright += '<p>Integration: Josh Thackray, 2015. <a href="mailto:thackray@rutgers.edu">thackray@rutgers.edu</a></p>'
    html =  '<!DOCTYPE html><html lang="en"><head><meta http-equiv="content-type" content="text/html; charset=UTF-8"><meta charset="utf-8"><title>'+filename+' .::. jvenn</title>'
    
    for ss in css:
        if ss.startswith('http'):
            html += '<link href="'+ss+'" rel="stylesheet" media="screen" />'
        else:
            html += '<style type="text/css">'+ss+'</style>'
            
    for js in js_src:
        html += '<script src="'+js+'"></script>'

    html += '<script language="Javascript">$(document).ready(function () { $("#jvenn-container").jvenn('+json.dumps(json_data)+'); });</script>'
    html += '</head><body><div class="container"><div class="row-fluid"><div class="span12"><div class="row-fluid"><div class="span7"><div class="row-fluid"><div id="jvenn-container"></div></div></div></div></div></div><hr><footer style="text-align: center; font-size:8px;">'
    html += cpyright+'</footer></div></body></html>'

    json_data = open(filename+'.html', "w")
    json_data.write(html)
    json_data.close()
    
    

if __name__ == "__main__":
    main()

