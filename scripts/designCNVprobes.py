import os
import sys
import time
import urllib2
import argparse
import requests
import subprocess
import pybedtools
import pybedtools.featurefuncs
from pybedtools import Interval
import xml.etree.ElementTree as ET
from ThackTech.isPCR import gfServer
from ThackTech.filetools import ensure_dir



#regions_filename = "test__finalvalidationSetFilteredSQNQ.bed" 
#regions_filename = "/mnt/data/projects/tourette/cnv_validation/finalvalidationSetFilteredSQNQ.bed"
#refseq_fasta_filename = ""
#subregions_per_region = 3
#primers_per_subregion = 10
#min_subregion_size = 1000
#genome_build = "hg19"
#genome_fasta = '/mnt/ref/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.2bit'
#cache_dir = 'data'

gargs = None


def get_arguments():
    global gargs
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #regions_filename
    parser.add_argument('bed', help="bed file describing the regions to design primers for")
    
    #genome_build
    parser.add_argument('--genome', default='hg19', help="UCSC genome build to use for designing primers. This should match the coordinate system used by your bed file")
    #genome_fasta
    parser.add_argument('--genome2bit', default='/mnt/ref/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.2bit', help="location of genome fasta in 2bit format")
    #subregions_per_region
    parser.add_argument('--subregions', type=int, default=3, help='Number of subregions to split each bed region into when looking for primers')
    #min_subregion_size
    parser.add_argument('--minsubregionsize', type=int, default=1000, help="Minimum size in bp for generated subregions")
    #primers_per_subregion
    parser.add_argument('--probespersubregion', type=int, default=10, help="Number of probes to design per subregion")
    
    
    #filter_group = parser.add_argument_group('Filtering Options')
    #filter_group.add_argument('')
    
    
    gfserver_group = parser.add_argument_group('isPCR server options')
    gfserver_group.add_argument('--gfserve', default="/home/josh/scripts/isPCR/gfServer", help="location of gfserver executable")
    gfserver_group.add_argument('--gfhost', default="localhost", help="hostname of the gfserver")
    gfserver_group.add_argument('--gfport', type=int, default=17779, help="port for the gfserver to use")
    
    #'hg38|17778|/mnt/ref/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.2bit'
    gfserver_group.add_argument('--altgfserve', action='append', default=[], help="Additional genomes to run isPCR on. specify as 'genome|port|reg-genome-2bit'")
    
    
    perf_group = parser.add_argument_group('Performance Options')
    #cache_dir
    perf_group.add_argument('--cachedir', default="data", help='directory to store cached files.')
    
    args = gargs = parser.parse_args()
    return args
#end get_arguments()

def main(args):
    gfservers = start_gfservers()
    ensure_dir(args.cachedir)
    regions = pybedtools.BedTool(args.bed)
    all_results = []
    filtered_all_results = []
    failed_all_results = []
    region_count = len(regions)
    for i in range(region_count):
        region = regions[i]
        print "Working on region {}/{} ({}; {}:{}-{}).....".format(i+1, region_count, region.name, region.chrom, region.start, region.stop)
        approx_subregion_size = region.length / args.subregions
        region_results = []
        if approx_subregion_size > args.minsubregionsize:
            print " -> splitting region into {} subregions of size {}....".format(args.subregions, approx_subregion_size)
            sub_regions = ["{}\t{}\t{}\t{}\n".format(region.chrom, 
                                                 region.start+(k*approx_subregion_size), 
                                                 region.start+((k+1)*approx_subregion_size), 
                                                 "{}.sr{}".format(region.name, k)
                                                ) for k in range(args.subregions)]
            sub_regions = pybedtools.BedTool("".join(sub_regions), from_string=True)
            
        else:
            print " -> region is too small ({}) for subregions. Processing entire region...".format(region.length)
            sub_regions = pybedtools.BedTool([region])
        
        for j in range(len(sub_regions)):
            subregion = sub_regions[j]
            region_results.extend(get_primers_for_region(subregion, args.probespersubregion))
        
        print " => Found {} total primer pairs for the entire region! <=".format(len(region_results))
        
        for j in range(len(region_results)):
            print " -> Validating result {} by in-silico PCR....".format(j+1)
            for server in gfservers:
                run_in_silico_pcr(region_results[j], server)
        
        #output the potential primers 
        make_region_candidate_bed(sub_regions, region_results, "{}_{}_{}-{}.primer_results.unfiltered.bed".format(region.name, region.chrom, region.start, region.stop))
        all_results.extend(region_results)
        
        #filter and output a more stringent set of primers
        filtered_region_results, failed_region_results = filter_candidates(region_results)
        make_region_candidate_bed(sub_regions, filtered_region_results, "{}_{}_{}-{}.primer_results.filtered.bed".format(region.name, region.chrom, region.start, region.stop))
        filtered_all_results.extend(filtered_region_results)
        failed_all_results.extend(failed_region_results)
        
        print "\n-----------------------------------------------------------\n"
    dump_results_to_file(sorted(filtered_all_results + failed_all_results, key=lambda x: x.name), "primer_results.unfiltered.tsv")
    dump_results_to_file(filtered_all_results, "primer_results.filtered.tsv")
    #stop_gfserver()    
#end main()



def filter_candidates(candidates):
    print " -> Filtering candidate primer pairs....."
    
    pass_results = []
    fail_results = []
    for pair in candidates:
        reasons = []
        
        for assembly in pair.ispcr:
            #Test that in-silico PCR found the intended amplicon
            if not pair.ispcr[assembly]['found_proper_amplicon']:
                reasons.append("       -> in-silico PCR could not find the intended amplicon in assembly {}.".format(assembly))
                #continue
                
            #Test that the off-target amplicon count is reasonable
            off_target_max = 0
            if pair.ispcr[assembly]['non_target_count'] > off_target_max:
                reasons.append("       -> in-silico PCR found {} off-target amplicons in assembly {}, more than max of {}.".format(pair.ispcr[assembly]['non_target_count'], assembly, off_target_max))
                #continue
        
            if len(pair.ispcr[assembly]['messages']) > 0:
                reasons.extend(pair.ispcr[assembly]['messages'])

        #check to make sure that amplicon covers the intended region
        if 'exon_overlap' not in pair.filter_exceptions:
            if len(pair.include_features['merged_intervals']) > 0 and not pybedtools.BedTool(pair.include_features['merged_intervals']).any_hits(pair.get_interval()):
                reasons.append("       -> Amplicon does not appear to overlap any exon features.")
                #continue
        
        #Test to make sure amplicon does not overlap known repetative regions from repeatmasker
        if fetch_repeats(pair.target_region).any_hits(pair.get_interval()):
            reasons.append("       -> Amplicon appears to overlap repeatitive elements.")
            #continue

            
        # Testing for GC clamp in primer sequence:
        # In the last 5 bases at the 3' end of the primer, make sure that there are at least 2 G or C bases (GC clamp).
        # G-C base pairs have a stronger bond than A-T base pairs (3 hydrogen bonds versus 2).
        min_num_gc_in_3prime_end = 2
        forward_3p_gc_count = pair.forward.sequence.count('C', -5) + pair.forward.sequence.count('G', -5)
        if forward_3p_gc_count < min_num_gc_in_3prime_end:
            reasons.append("       -> Forward primer contains too few ({}) G/C bases in the 3' end.".format(forward_3p_gc_count))
            
        reverse_3p_gc_count = pair.reverse.sequence.count('C', -5) + pair.reverse.sequence.count('G', -5)
        if reverse_3p_gc_count < min_num_gc_in_3prime_end:
            reasons.append("       -> Reverse primer contains too few ({}) G/C bases in the 3' end.".format(reverse_3p_gc_count))
        
            
            
        #filter on total number of SNPs present in the priming region (not including the amplified region)
        region_snps = fetch_snps(pair.target_region)
        
        max_snps_in_primer = 2
        forward_snp_hits = region_snps.all_hits(pair.forward.get_interval())
        reverse_snp_hits = region_snps.all_hits(pair.reverse.get_interval())
        if len(forward_snp_hits) > max_snps_in_primer:
            reasons.append("       -> Forward primer contains too many ({}) SNPs ({}).".format(len(forward_snp_hits), ";".join([hit.name for hit in forward_snp_hits])))
            
        if len(reverse_snp_hits) > max_snps_in_primer:
            reasons.append("       -> Reverse primer contains too many ({}) SNPs ({}).".format(len(reverse_snp_hits), ";".join([hit.name for hit in reverse_snp_hits])))
        
        
        #filter on any SNP presence in the 3' end of the primer
        max_snps_in_3p_primer = 0
        num_3p_bases_for_snp_check = 5
        forward_3p_snp_hits = region_snps.all_hits(pybedtools.featurefuncs.three_prime(pair.forward.get_interval(), upstream=num_3p_bases_for_snp_check, downstream=0))
        reverse_3p_snp_hits = region_snps.all_hits(pybedtools.featurefuncs.three_prime(pair.reverse.get_interval(), upstream=num_3p_bases_for_snp_check, downstream=0))
        if len(forward_3p_snp_hits) > max_snps_in_3p_primer:
            reasons.append("       -> Forward primer contains {} SNPs in 3' region ({}).".format(len(forward_3p_snp_hits), ";".join([hit.name for hit in forward_3p_snp_hits])))
            
        if len(reverse_3p_snp_hits) > max_snps_in_3p_primer:
            reasons.append("       -> Reverse primer contains {} SNPs in 3' region ({}).".format(len(reverse_3p_snp_hits), ";".join([hit.name for hit in reverse_3p_snp_hits])))
        
        pair.filter_results = reasons
        if len(reasons) > 0:
            print "    -> Excluding pair {} for the following reasons:".format(pair.name)
            for reason in reasons:
                print reason
            fail_results.append(pair)
        else:
            #passed all filtering criteria
            pass_results.append(pair)

    print "    => Filtering results: {} (pass); {} (fail); {} (total examined)".format(len(pass_results), len(fail_results), len(candidates))
    return (pass_results, fail_results)
#end filter_candidates()



def make_region_candidate_bed(regions, results, filename, include_headers=True):
    pos_min = sys.maxint
    pos_max = 0
    for r in regions:
        if r.start < pos_min:
            pos_min = r.start
        if r.stop > pos_max:
            pos_max = r.stop

    with open(filename, 'w+') as bed_file:
        bed_file.write('browser position {}:{}-{}\n'.format(regions[0].chrom, pos_min, pos_max))
        bed_file.write('track db="{}" name="{}" visibility=2\n'.format(gargs.genome, "Region Canidate Primers"))
        bed_file.write("{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t255,0,0\t1\t{}\t0\n".format(regions[0].chrom, pos_min, pos_max, regions[0].name.replace('.sr0', ""), pos_min, pos_max, pos_max - pos_min))
        for region in regions:
            bed_file.write("{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t0,255,0\t1\t{}\t0\n".format(region.chrom, region.start, region.stop, region.name, region.start, region.stop, region.length))
        for result in results:
            bed_file.write(result.get_bed_line())
            bed_file.write("\n")
#end make_region_candidate_bed()


def dump_results_to_file(results, filename):
    
    with open(filename, 'w+') as f:
        f.write(("{}\t{}\t{}\t"
                 + "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t"
				 + "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t")
               .format("Pair_Name", "Product_Size", "Overall_Penalty",
                   "Forward_Name", "Forward_Sequence", "Forward_Start", "Forward_Stop", "Forward_Length", "Forward_Penalty", "Forward_Tm", "Forward_GC_Percent",
				   "Reverse_Name", "Reverse_Sequence", "Reverse_Start", "Reverse_Stop", "Reverse_Length", "Reverse_Penalty", "Reverse_Tm", "Reverse_GC_Percent"
                   )
                )
        for assembly in results[0].ispcr:
            f.write("{}\t{}\t{}\t".format("isPCR["+assembly+"]_amplicon_count", "isPCR["+assembly+"]_found_target", "isPCR["+assembly+"]_off_target_count"))
        
        f.write("{}\t".format("Filtering_Results"))
        f.write("\n")
        for pair in results:
            f.write("{}\t{}\t{}\t".format(pair.name, pair.size, pair.penalty))
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(pair.forward.name, pair.forward.sequence.upper(), pair.forward.start, pair.forward.stop, pair.forward.length, pair.forward.penalty, pair.forward.tm, pair.forward.gc))
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(pair.reverse.name, pair.reverse.sequence.upper(), pair.reverse.start, pair.reverse.stop, pair.reverse.length, pair.reverse.penalty, pair.reverse.tm, pair.reverse.gc))
            for assembly in pair.ispcr:
                f.write('{}\t{}\t{}\t'.format(pair.ispcr[assembly]['count'], pair.ispcr[assembly]['found_proper_amplicon'], pair.ispcr[assembly]['non_target_count']))
            
            f.write('{}\t'.format("; ".join(pair.filter_results)))
            f.write('\n')
#end dump_results_to_file()
            
def fetch_sequence(region):
    cache_path = ".{}.sequence.segment={}.{}-{}.fasta".format(gargs.genome, region.chrom, region.start, region.stop)
    cache_path = os.path.join(gargs.cachedir, cache_path)
    if not os.path.exists(cache_path):
        f = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/{}/dna?segment={}:{},{}".format(gargs.genome, region.chrom, region.start, region.stop))
        root = ET.fromstring(f.read())
        f.close()
        with open(cache_path, 'w+') as cf:
            cf.write(root.findtext("./SEQUENCE/DNA").replace('\n','').replace('\r',''))
    
    with open(cache_path, 'r') as f:
        data = f.read()
    return data
#end fetch_sequence()

def fetch_repeats(region):
    cache_path = ".{}.repeats.{}.{}-{}.bed".format(gargs.genome, region.chrom, region.start, region.stop)
    cache_path = os.path.join(gargs.cachedir, cache_path)
    if not os.path.exists(cache_path):
        cmd = [
            'mysql',
            '--user', 'genome',
            '--host', 'genome-mysql.cse.ucsc.edu',
            '-N', '-AB', 
            '-e', 'SELECT genoName, genoStart, genoEnd, concat(repClass, "_", repFamily, "_", repName) as name FROM {}.rmsk WHERE genoName = "{}" and (genoStart > {} or genoEnd < {});'.format(gargs.genome, region.chrom, region.start, region.stop),
            gargs.genome
        ]
        #print " ".join(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdoutdata, stderrdata) = p.communicate()
        repeats = pybedtools.BedTool(stdoutdata, from_string=True)
        repeats.saveas(cache_path)
    else:
        repeats = pybedtools.BedTool(cache_path)
    
    return repeats;
#end fetch_repeats()

def fetch_features(region):
    print "     -> Fetching gene features..."
    #server = "http://grch37.rest.ensembl.org" #provides GRCh37 data
    #server = "http://rest.ensembl.org" #only provides GRCh38 build data!!!!!
    cache_path = ".{}.features.segment={}.{}-{}.type={}".format(gargs.genome, region.chrom, region.start, region.stop, "refGene")
    cache_path = os.path.join(gargs.cachedir, cache_path)
    if not os.path.exists(cache_path):
        f = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/{}/features?segment={}:{},{};type={}".format(gargs.genome, region.chrom, region.start, region.stop, "refGene")) #"ccdsGene"))
        with open(cache_path, 'w+') as cf:
            cf.write(f.read())
        f.close()
        
    with open(cache_path, 'r') as f:
        root = ET.fromstring(f.read())
        
    
    all_features = root.findall("./GFF/SEGMENT/FEATURE")
    unique_labels = set([feature.get("label") for feature in all_features])
    
    reportable_features = []
    for label in unique_labels:
        group = root.findall("./GFF/SEGMENT/FEATURE[@label='{}']".format(label))
        record = {
            'id': label,
            'orientation': group[0].findtext('ORIENTATION'),
            'exons': [],
        }
        for item in group:
            record['exons'].append(Interval(region.chrom, int(item.findtext('START')), int(item.findtext('END')), item.get('id'), strand=item.findtext('ORIENTATION')))
        reportable_features.append(record)
    
    if len(reportable_features) > 0:
        merged = pybedtools.BedTool([item for sublist in reportable_features for item in sublist['exons']]).sort().merge().intersect(pybedtools.BedTool([region]), u=True)
        print "       -> Found {} gene features containg {} exonic regions".format(len(reportable_features), len(merged))
        if len(merged) > 0:
            complement = merged.complement(genome=gargs.genome).intersect(pybedtools.BedTool([region]), u=True).saveas()#.filter(lambda b: b.chrom == region.chrom and b.start > region.start and b.stop < region.stop).saveas()
        else:
            complement = [region]
    else:
        print "       -> Found 0 gene features"
        merged = []
        complement = [region]
    
    return {
        'count': len(unique_labels),
        'features': reportable_features,
        'merged_intervals': list(merged),
        'complement_intervals': list(complement) 
    }
#end fetch_features()


def fetch_snps(region):
    cache_path = ".{}.snps.{}.{}-{}.bed".format(gargs.genome, region.chrom, region.start, region.stop)
    cache_path = os.path.join(gargs.cachedir, cache_path)
    if not os.path.exists(cache_path):
        cmd = [
            'mysql',
            '--user', 'genome',
            '--host', 'genome-mysql.cse.ucsc.edu',
            '-N', '-AB', 
            '-e', 'SELECT chrom, chromStart, chromEnd, CONCAT(name, "[", class, "][", observed, "]")  '
                + 'FROM {}.snp147Common '
                + 'WHERE chrom = "{}" AND ((chromStart BETWEEN {} AND {}) OR (chromEnd BETWEEN {} AND {})) AND molType = "genomic";'.format(gargs.genome, region.chrom, region.start, region.stop, region.start, region.stop),
            gargs.genome
        ]
        #print " ".join(cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdoutdata, stderrdata) = p.communicate()
        snps = pybedtools.BedTool(stdoutdata, from_string=True)
        snps.saveas(cache_path)
    else:
        snps = pybedtools.BedTool(cache_path)
    
    return snps;
#end fetch_snps()


def get_primers_for_region(region, num_primers=5, override_settings={}):
    ensure_dir('primer3')
    print "   -> Looking for primers in region {}:{}-{}".format(region.chrom, region.start, region.stop)
    settings = get_default_primer3_settings()
    settings['SEQUENCE_ID'] = "{}_{}_{}-{}".format(region.name, region.chrom, region.start, region.stop)
    settings['SEQUENCE_TEMPLATE'] = fetch_sequence(region)
    settings['PRIMER_NUM_RETURN'] = num_primers
    #settings['SEQUENCE_INCLUDED_REGION']
    #settings['SEQUENCE_TARGET'] 
    settings['SEQUENCE_EXCLUDED_REGION'] = []#'1,{}'.format(region.length-1,)
    
    filter_exceptions = []
     
    include_features = fetch_features(region)
    if len(include_features['merged_intervals']) > 0:
        for f in include_features['complement_intervals']:
            s = f.start - region.start
            e = f.stop - region.start
            #l = f.length
            if s < 1:
                #l += s
                s = 1
            if e >= region.length:
                e = region.length - 1
            settings['SEQUENCE_EXCLUDED_REGION'].append("{},{}".format(s, e-s))
        #settings['SEQUENCE_TARGET'] = ",".join(settings['SEQUENCE_TARGET'])
    else:
        filter_exceptions.append('exon_overlap')
        print "         -> It looks like there are no gene features in this region, so not including exon overlap constraint in primer picking."
    
    settings.update(override_settings)
    
    results, errors, warnings = run_primer3(region, settings)
    
    
    if len(results) == 0:  #implement relaxation of settings and retry
        if errors == "Too many elements for tag SEQUENCE_EXCLUDED_REGION":
            print "         -> Attempting to relax parameters for search:"
            print "             -> Noticed that there are **MANY** exons in this region; removing exon overlap constraint (will be checked during filtering) and resubmitting (see next entry --v)"
            settings.update({'SEQUENCE_EXCLUDED_REGION': []})
            results, errors, warnings = run_primer3(region, settings)
        
        elif len(include_features['merged_intervals']) < 3:
            print "         -> Attempting to relax parameters for search:"
            print "             -> Noticed that there <3 exon(s) in this region; removing exon overlap constraint and resubmitting (see next entry --v)"
            settings.update({'SEQUENCE_EXCLUDED_REGION': []})
            filter_exceptions.append('exon_overlap')
            results, errors, warnings = run_primer3(region, settings)
        
    
    for r in results:
        r.filter_exceptions = filter_exceptions
        r.target_region = region
        r.include_features = include_features
    
    return results
#end get_primers_for_region()

def run_primer3(region, settings):
    p = subprocess.Popen(['primer3_core'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    indata = ""
    for (key,val) in settings.iteritems():
        if isinstance(val, list):
            indata += "{}={}\n".format(key, " ".join(val))
        else:
            indata += "{}={}\n".format(key,val)
    indata += "=\n"
    with open("primer3/primer3.{}_{}_{}-{}.params.txt".format(region.name, region.chrom, region.start, region.stop), 'w+') as cf:
        cf.write(indata)
    #print indata
    (stdoutdata, stderrdata) = p.communicate(indata)
    with open("primer3/primer3.{}_{}_{}-{}.output.txt".format(region.name, region.chrom, region.start, region.stop), 'w+') as cf:
        cf.write(stdoutdata)
        cf.write(stderrdata)
    #print stdoutdata
    #print stderrdata
    results, errors, warnings = parse_primer3_output(region, stdoutdata)
    
    return results, errors, warnings
#end run_primer3()
    

def get_default_primer3_settings():
    return {
        'SEQUENCE_ID': '',
        'SEQUENCE_TEMPLATE': '',
#        'SEQUENCE_TARGET':    '',
        'PRIMER_LOWERCASE_MASKING': 0,

        'PRIMER_MIN_THREE_PRIME_DISTANCE': 100, #force primers returned to be at least 100bp apart
        
        'PRIMER_NUM_RETURN': 5, #number of primers to return
        
        'PRIMER_TASK': 'generic',            #We want just normal PCR primers
        'PRIMER_PICK_LEFT_PRIMER': 1,        #Pick a Left oligo
        'PRIMER_PICK_INTERNAL_OLIGO': 0,    #do not pick an internal oligo
        'PRIMER_PICK_RIGHT_PRIMER': 1,        #pick a right oligo
        
        'PRIMER_PRODUCT_SIZE_RANGE': '80-180',    #return primers that generate amplicons in this range
        'PRIMER_OPT_SIZE': 20,    #primer optimum length
        'PRIMER_MIN_SIZE': 18,    #primer minimum length
        'PRIMER_MAX_SIZE': 27,    #primer maximum length
        
        'PRIMER_OPT_TM': 60,    #primer optimum melting temp
        'PRIMER_MIN_TM': 57,    #primer minimum melting temp
        'PRIMER_MAX_TM': 63,    #primer maximum melting temp
        'PRIMER_PAIR_MAX_DIFF_TM': 2,    #melting temperature should not differ by more than 2 degrees
        
        'PRIMER_OPT_GC_PERCENT': 50,    #primer optimum GC%
        'PRIMER_MIN_GC': 40,            #primer minimum GC%
        'PRIMER_MAX_GC': 60,            #primer maximum GC%
        
        'PRIMER_MAX_SELF_ANY': 2,
        'PRIMER_MAX_SELF_END': 0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 4,
        'PRIMER_PAIR_MAX_COMPL_END': 0,
        
        'PRIMER_MISPRIMING_LIBRARY': '/home/josh/scripts/isPCR/misprime_lib/human_rep_and_simple.fa', #mispriming library to screen against
        'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 1,
        'PRIMER_MAX_LIBRARY_MISPRIMING': 12,
        'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 24,
        
        'PRIMER_MAX_NS_ACCEPTED': '1', #Do not allow primer to contain more than 1 ambiguous base
        'PRIMER_MAX_POLY_X': 3,    #Allow a maximum of tri-nucleotide repeat
        
#        'P3_FILE_FLAG':    '1',
#        'PRIMER_EXPLAIN_FLAG':    '1'        
    }
#end get_default_primer3_settings()

def parse_primer3_output(region, data):
    #print data
    data_dict = dict(line.split("=") for line in data.split("\n") if line  not in ["=", ""])
    errors = None
    warnings = None
    if 'PRIMER_ERROR' in data_dict:
        errors = data_dict['PRIMER_ERROR']
        print "         => Primer3 ERROR: {} <=".format(data_dict['PRIMER_ERROR'])
        print "         => Returning no primers for this region! <="
        return [], errors, warnings
        
    if 'PRIMER_WARNING' in data_dict:
        warnings = data_dict['PRIMER_WARNING']
        print "         => Primer3 WARNING: {} <=".format(data_dict['PRIMER_WARNING'])
    
    num_results = int(data_dict['PRIMER_PAIR_NUM_RETURNED'])
    results = []
    if num_results > 0:
        for i in range(num_results):
            forward = Primer(gargs.genome,
							 "LEFT", 
                             "{}_{}_{}".format(data_dict['SEQUENCE_ID'], i, "F"),
                             data_dict['PRIMER_LEFT_{}_SEQUENCE'.format(i)],
                             region.chrom,
                             region.start + int(data_dict['PRIMER_LEFT_{}'.format(i)].split(',')[0]) - 1,
                             int(data_dict['PRIMER_LEFT_{}'.format(i)].split(',')[1]),
                             '+', # "LEFT" primer is always 5' -> 3' on the same strand as the input sequence
                             float(data_dict['PRIMER_LEFT_{}_PENALTY'.format(i)]),
                             float(data_dict['PRIMER_LEFT_{}_TM'.format(i)]),
                             float(data_dict['PRIMER_LEFT_{}_GC_PERCENT'.format(i)]),
                             dict((key.replace('PRIMER_LEFT_{}'.format(i), ''),val) for (key,val) in data_dict.iteritems() if key.startswith('PRIMER_LEFT_{}'.format(i)))
                            )
            reverse = Primer(gargs.genome,
							 "RIGHT", 
                             "{}_{}_{}".format(data_dict['SEQUENCE_ID'], i, "R"),
                             data_dict['PRIMER_RIGHT_{}_SEQUENCE'.format(i)],
                             region.chrom,
                             region.start + int(data_dict['PRIMER_RIGHT_{}'.format(i)].split(',')[0]) - int(data_dict['PRIMER_RIGHT_{}'.format(i)].split(',')[1]),
                             int(data_dict['PRIMER_RIGHT_{}'.format(i)].split(',')[1]),
                             '-', #  "RIGHT" primer is presented 5' -> 3' on the opposite strand from the input sequence
                             float(data_dict['PRIMER_RIGHT_{}_PENALTY'.format(i)]),
                             float(data_dict['PRIMER_RIGHT_{}_TM'.format(i)]),
                             float(data_dict['PRIMER_RIGHT_{}_GC_PERCENT'.format(i)]),
                             dict((key.replace('PRIMER_RIGHT_{}'.format(i), ''),val) for (key,val) in data_dict.iteritems() if key.startswith('PRIMER_RIGHT_{}'.format(i)))
                            )
            pair = PrimerResult(gargs.genome,
								"{}_{}".format(data_dict['SEQUENCE_ID'], i), 
                                forward, 
                                reverse, 
                                int(data_dict['PRIMER_PAIR_{}_PRODUCT_SIZE'.format(i)]),
                                float(data_dict['PRIMER_PAIR_{}_PENALTY'.format(i)]),
                                dict((key.replace('PRIMER_PAIR_{}'.format(i), ''),val) for (key,val) in data_dict.iteritems() if key.startswith('PRIMER_PAIR_{}'.format(i)))
                                )
            results.append(pair)
        print "         -> Found {} primers for this region!".format(num_results)
    else:
        print "         => No primers found for this region! <="
    return (results, errors, warnings)
#end parse_primer3_output()


def start_gfservers():
    servers = []
    servers.append(gfServer(gargs.genome, gargs.gfserve, gargs.gfhost, gargs.gfport, gargs.genome2bit))
    
    for i in range(len(gargs.altgfserve)):
        parts = gargs.altgfserve[i].split('|')
        servers.append(gfServer(parts[0], gargs.gfserve, gargs.gfhost, int(parts[1]), parts[2]))
    
    for server in servers:
        server.start()
        
    return servers
#end start_gfservers()

def run_in_silico_pcr(primerpair, server):
    ensure_dir('ispcr')
    hits = server.isPCR(primerpair.forward.sequence, primerpair.reverse.sequence, out='bed', maxsize=4000, minPerfect=15, minGood=15)
    hits.saveas("ispcr/{}.ispcr.{}.bed".format(primerpair.name, server.name))
    
    count = len(hits)
    pair_interval = primerpair.get_interval(server.name)
    messages = []
    if pair_interval is None:
        print "Did not get an interval, so cannot assess specificity of isPCR results (will assume all off-target), but found {} total amplicons...".format(count)
        found_proper_amplicon = False
        non_target_count = count
        messages.append('Liftover of interval from assembly {} to assembly {} failed, so cannot assess specificity of isPCR results (will assume all off-target)'.format(primerpair.assembly, server.name))
    else:
        found_proper_amplicon = hits.any_hits(pair_interval)
        non_target_count = len(hits) - int(hits.any_hits(pair_interval))
        print "   -> Found {} total amplicons, {}including targeted site, and {} off-target amplicons in assembly {}".format(count, ("" if found_proper_amplicon else "NOT "), non_target_count, server.name)

    if not hasattr(primerpair, 'ispcr'):
        primerpair.ispcr = {}
        
    primerpair.ispcr[server.name] = {
		'assembly': server.name,
        'count': count,
        'found_proper_amplicon': found_proper_amplicon,
        'non_target_count': non_target_count,
		'hits': hits,
		'messages': messages
    }

#end run_in_silico_pcr()


###############################
#                             #
# Classes and Data Structures #
#                             #
###############################
class PrimerResult:

    def __init__(self, assembly, name, forward, reverse, size, penalty, props):
        self.assembly = assembly
        self.name = name
        self.forward = forward
        self.reverse = reverse
        self.size = size
        self.penalty = penalty
        self.properties = props
    #end __init__()
    
    def get_bed_line(self, color=(0,0,255)):
        return "{}\t{}\t{}\t{}\t0\t+\t{}\t{}\t{}\t{}\t{}\t{}".format(self.forward.chrom, 
                                                                 self.forward.start, 
                                                                 self.reverse.stop, 
                                                                 self.name, 
                                                                 self.forward.start, 
                                                                 self.reverse.stop, 
                                                                 ",".join(color), 
                                                                 2, 
                                                                 ",".join([str(self.forward.length), str(self.reverse.length)]), 
                                                                 ",".join([str(self.forward.start - self.forward.start), str(self.reverse.start - self.forward.start)]))
    def get_interval(self, assembly=None):
        i = Interval(self.forward.chrom, min(self.forward.start, self.reverse.start), max(self.forward.stop, self.reverse.stop), self.name, strand='+')
        i.file_type = 'bed'
        if assembly is not None and assembly != self.assembly:
            i = liftover(self.assembly, i, assembly)
        return i
#end class PrimerResult

class Primer:
    def __init__(self, assembly, orientation, name, seq, chrom, start, length, strand, penalty, tm, gc, props):
        self.assembly = assembly
        self.orientation = orientation
        self.name = name
        self.sequence = seq.upper()
        self.chrom = chrom
        self.start = start
        self.length = length
        self.stop = self.start + self.length
        self.strand = strand
        self.penalty = penalty
        self.tm = tm
        self.gc = gc
        self.properties = dict(props)

    def get_bed_line(self):
        return "{}\t{}\t{}\t{}\t0\t+".format(self.chrom, self.start, self.stop, self.name)
    
    def get_interval(self, assembly=None):
        i = Interval(self.chrom, self.start, self.stop, self.name, strand=self.strand)
        i.file_type = 'bed'
        if assembly is not None and assembly != self.assembly:
            i = liftover(self.assembly, i, assembly)
        return i
#end class Primer

def liftover(curr_assembly, interval, target_assembly, species='homo_sapiens'):
    ass_names = {
        'hg19': 'GRCh37',
        'hg38': 'GRCh38'
    }
    if curr_assembly in ass_names:
        curr_assembly = ass_names[curr_assembly]
    if target_assembly in ass_names:
        target_assembly = ass_names[target_assembly]    

    server = "http://grch37.rest.ensembl.org"
    ext = "/map/{}/{}/{}:{}..{}:{}/{}?".format(species, curr_assembly, interval.chrom.replace('chr', ''), interval.start, interval.stop, -1 if interval.strand == '-' else 1, target_assembly)
    #print server+ext
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
        r.raise_for_status()
        sys.exit()
     
    decoded = r.json()
    #print repr(decoded)
    if len(decoded['mappings']) > 0:
        lifted = decoded['mappings'][0]['mapped']
        #print lifted
        liftedInterval = Interval('chr'+lifted['seq_region_name'], lifted['start'], lifted['end'], interval.name+'_'+target_assembly, strand=('-' if lifted['strand'] == -1 else '+'))
        liftedInterval.file_type = 'bed'
        return liftedInterval
    else:
        print "Unable to convert interval {}:{}-{} from assembly {} to assembly {}".format(interval.chrom, interval.start, interval.stop, curr_assembly, target_assembly)
        return None
#end liftover()

###################################
#                                 #
# END Classes and Data Structures #
#                                 #
###################################





if __name__ == "__main__":
    main(get_arguments())
        
        
        
        



#########################################
#
# This is the entire list of parameters  
# that Primer3 accepts
#
#########################################
# SEQUENCE_EXCLUDED_REGION
# SEQUENCE_INCLUDED_REGION
# SEQUENCE_PRIMER_REVCOMP
# SEQUENCE_FORCE_LEFT_END
# SEQUENCE_INTERNAL_EXCLUDED_REGION
# SEQUENCE_QUALITY
# SEQUENCE_FORCE_LEFT_START
# SEQUENCE_INTERNAL_OLIGO
# SEQUENCE_START_CODON_POSITION
# SEQUENCE_FORCE_RIGHT_END
# SEQUENCE_OVERLAP_JUNCTION_LIST
# SEQUENCE_TARGET
# SEQUENCE_FORCE_RIGHT_START
# SEQUENCE_PRIMER
# SEQUENCE_TEMPLATE
# SEQUENCE_ID
# SEQUENCE_PRIMER_PAIR_OK_REGION_LIST

# PRIMER_DNA_CONC
# PRIMER_MAX_END_GC
# PRIMER_PAIR_WT_PRODUCT_SIZE_LT
# PRIMER_DNTP_CONC
# PRIMER_MAX_END_STABILITY
# PRIMER_PAIR_WT_PRODUCT_TM_GT
# PRIMER_EXPLAIN_FLAG
# PRIMER_MAX_GC
# PRIMER_PAIR_WT_PRODUCT_TM_LT
# PRIMER_FIRST_BASE_INDEX
# PRIMER_MAX_HAIRPIN_TH
# PRIMER_PAIR_WT_PR_PENALTY
# PRIMER_GC_CLAMP
# PRIMER_MAX_LIBRARY_MISPRIMING
# PRIMER_PAIR_WT_TEMPLATE_MISPRIMING
# PRIMER_INSIDE_PENALTY
# PRIMER_MAX_NS_ACCEPTED
# PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH
# PRIMER_INTERNAL_DNA_CONC
# PRIMER_MAX_POLY_X
# PRIMER_PICK_ANYWAY
# PRIMER_INTERNAL_DNTP_CONC
# PRIMER_MAX_SELF_ANY
# PRIMER_PICK_INTERNAL_OLIGO
# PRIMER_INTERNAL_MAX_GC
# PRIMER_MAX_SELF_ANY_TH
# PRIMER_PICK_LEFT_PRIMER
# PRIMER_INTERNAL_MAX_HAIRPIN_TH
# PRIMER_MAX_SELF_END
# PRIMER_PICK_RIGHT_PRIMER
# PRIMER_INTERNAL_MAX_LIBRARY_MISHYB
# PRIMER_MAX_SELF_END_TH
# PRIMER_PRODUCT_MAX_TM
# PRIMER_INTERNAL_MAX_NS_ACCEPTED
# PRIMER_MAX_SIZE
# PRIMER_PRODUCT_MIN_TM
# PRIMER_INTERNAL_MAX_POLY_X
# PRIMER_MAX_TEMPLATE_MISPRIMING
# PRIMER_PRODUCT_OPT_SIZE
# PRIMER_INTERNAL_MAX_SELF_ANY
# PRIMER_MAX_TEMPLATE_MISPRIMING_TH
# PRIMER_PRODUCT_OPT_TM
# PRIMER_INTERNAL_MAX_SELF_ANY_TH
# PRIMER_MAX_TM
# PRIMER_PRODUCT_SIZE_RANGE
# PRIMER_INTERNAL_MAX_SELF_END
# PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION
# PRIMER_QUALITY_RANGE_MAX
# PRIMER_INTERNAL_MAX_SELF_END_TH
# PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION
# PRIMER_QUALITY_RANGE_MIN
# PRIMER_INTERNAL_MAX_SIZE
# PRIMER_MIN_END_QUALITY
# PRIMER_SALT_CORRECTIONS
# PRIMER_INTERNAL_MAX_TM
# PRIMER_MIN_GC
# PRIMER_SALT_DIVALENT
# PRIMER_INTERNAL_MIN_GC
# PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE
# PRIMER_SALT_MONOVALENT
# PRIMER_INTERNAL_MIN_QUALITY
# PRIMER_MIN_QUALITY
# PRIMER_SEQUENCING_ACCURACY
# PRIMER_INTERNAL_MIN_SIZE
# PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE
# PRIMER_SEQUENCING_INTERVAL
# PRIMER_INTERNAL_MIN_TM
# PRIMER_MIN_SIZE
# PRIMER_SEQUENCING_LEAD
# PRIMER_INTERNAL_MISHYB_LIBRARY
# PRIMER_MIN_THREE_PRIME_DISTANCE
# PRIMER_SEQUENCING_SPACING
# PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME
# PRIMER_MIN_TM
# PRIMER_TASK
# PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME
# PRIMER_MISPRIMING_LIBRARY
# PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT
# PRIMER_INTERNAL_OPT_GC_PERCENT
# PRIMER_MUST_MATCH_FIVE_PRIME
# PRIMER_THERMODYNAMIC_PARAMETERS_PATH
# PRIMER_INTERNAL_OPT_SIZE
# PRIMER_MUST_MATCH_THREE_PRIME
# PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT
# PRIMER_INTERNAL_OPT_TM
# PRIMER_NUM_RETURN
# PRIMER_TM_FORMULA
# PRIMER_INTERNAL_SALT_DIVALENT
# PRIMER_OPT_GC_PERCENT
# PRIMER_WT_END_QUAL
# PRIMER_INTERNAL_SALT_MONOVALENT
# PRIMER_OPT_SIZE
# PRIMER_WT_END_STABILITY
# PRIMER_INTERNAL_WT_END_QUAL
# PRIMER_OPT_TM
# PRIMER_WT_GC_PERCENT_GT
# PRIMER_INTERNAL_WT_GC_PERCENT_GT
# PRIMER_OUTSIDE_PENALTY
# PRIMER_WT_GC_PERCENT_LT
# PRIMER_INTERNAL_WT_GC_PERCENT_LT
# PRIMER_PAIR_MAX_COMPL_ANY
# PRIMER_WT_HAIRPIN_TH
# PRIMER_INTERNAL_WT_HAIRPIN_TH
# PRIMER_PAIR_MAX_COMPL_ANY_TH
# PRIMER_WT_LIBRARY_MISPRIMING
# PRIMER_INTERNAL_WT_LIBRARY_MISHYB
# PRIMER_PAIR_MAX_COMPL_END
# PRIMER_WT_NUM_NS
# PRIMER_INTERNAL_WT_NUM_NS
# PRIMER_PAIR_MAX_COMPL_END_TH
# PRIMER_WT_POS_PENALTY
# PRIMER_INTERNAL_WT_SELF_ANY
# PRIMER_PAIR_MAX_DIFF_TM
# PRIMER_WT_SELF_ANY
# PRIMER_INTERNAL_WT_SELF_ANY_TH
# PRIMER_PAIR_MAX_LIBRARY_MISPRIMING
# PRIMER_WT_SELF_ANY_TH
# PRIMER_INTERNAL_WT_SELF_END
# PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
# PRIMER_WT_SELF_END
# PRIMER_INTERNAL_WT_SELF_END_TH
# PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH
# PRIMER_WT_SELF_END_TH
# PRIMER_INTERNAL_WT_SEQ_QUAL
# PRIMER_PAIR_WT_COMPL_ANY
# PRIMER_WT_SEQ_QUAL
# PRIMER_INTERNAL_WT_SIZE_GT
# PRIMER_PAIR_WT_COMPL_ANY_TH
# PRIMER_WT_SIZE_GT
# PRIMER_INTERNAL_WT_SIZE_LT
# PRIMER_PAIR_WT_COMPL_END
# PRIMER_WT_SIZE_LT
# PRIMER_INTERNAL_WT_TM_GT
# PRIMER_PAIR_WT_COMPL_END_TH
# PRIMER_WT_TEMPLATE_MISPRIMING
# PRIMER_INTERNAL_WT_TM_LT
# PRIMER_PAIR_WT_DIFF_TM
# PRIMER_WT_TEMPLATE_MISPRIMING_TH
# PRIMER_LIBERAL_BASE
# PRIMER_PAIR_WT_IO_PENALTY
# PRIMER_WT_TM_GT
# PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
# PRIMER_PAIR_WT_LIBRARY_MISPRIMING
# PRIMER_WT_TM_LT
# PRIMER_LOWERCASE_MASKING
# PRIMER_PAIR_WT_PRODUCT_SIZE_GT