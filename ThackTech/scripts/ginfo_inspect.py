#!/usr/bin/env python

import os
import sys
import glob
import argparse
import tabulate
from ThackTech import conf
from ThackTech.Pipelines import GenomeInfo



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('genome', nargs='?', help="Name of reference genome to inspect.")
	parser.add_argument('-c', '--config', help="specify additional genome configuration file.")
	args = parser.parse_args()
	
	if args.config is not None:
		sys.stderr.write('Registering configuration at "{}"\n\n'.format(args.config))
		conf.register_config_location('genomes', args.config)
	
	
	ref_genomes = GenomeInfo.get_reference_genomes()
	
	if args.genome is not None and args.genome not in ref_genomes:
		sys.stderr.write('Could not find any entry for genome {}!\n'.format(args.genome))
		sys.exit(1)
	
	for gname in ref_genomes.keys():
		g = ref_genomes[gname]
		
		if args.genome is not None and args.genome != gname:
			continue #continue if not specific genome

		sys.stderr.write('Inspecting Genome "{}":\n'.format(gname))
		
		data = []
		data.append({
			'Attribute': 'Genome Size',
			'Value': g.gsize,
			'Exists?': 'N/A'
		})
		data.append({
			'Attribute': 'Chrom Sizes',
			'Value': g.chrsize,
			'Exists?': bool_to_str(file_exists(g.chrsize))
		})
		data.append({
			'Attribute': 'Whole Genome FASTA',
			'Value': g.wg_fasta,
			'Exists?': bool_to_str(file_exists(g.wg_fasta))
		})
		for chrom_fasta in sorted(g.chr_fasta):
			data.append({
				'Attribute': '{} FASTA'.format(chrom_fasta),
				'Value': g.chr_fasta[chrom_fasta],
				'Exists?': bool_to_str(file_exists(g.chr_fasta[chrom_fasta]))
			})
		data.append({
			'Attribute': 'Gene Annotation File',
			'Value': g.genes_gtf,
			'Exists?': bool_to_str(file_exists(g.genes_gtf))
		})
		for idx in sorted(g.indexes):
			data.append({
				'Attribute': idx,
				'Value': g.indexes[idx],
				'Exists?': str(len(glob.glob(g.indexes[idx]+'*')) > 0)
			})
		
		sys.stderr.write(tabulate.tabulate(data, headers='keys', tablefmt="grid"))
		sys.stderr.write('\n\n\n')
#end main()
	
def bool_to_str(b):
	if b:
		return 'Yes'
	else:
		return 'No'
#end bool_to_str()

def file_exists(path):
	if path is None:
		return False
	if os.path.exists(path) and os.path.isfile(path):
		return True
	return False
#end file_exists()



if __name__ == "__main__":
	main()










	
