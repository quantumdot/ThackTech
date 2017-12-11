#!/usr/bin/env python

import os
import sys
import glob
import argparse
import tabulate
from ThackTech.Pipelines import GenomeInfo



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('genome', nargs='?', help="Name of reference genome to inspect.")
	
	args = parser.parse_args()
	
	
	ref_genomes = GenomeInfo.get_reference_genomes()
	
	if args.genome is not None and args.genome not in ref_genomes:
		sys.stderr.writ('Could not find any entry for genome {}!\n'.format(args.genome))
		sys.exit(1)
	
	for gname in ref_genomes.keys():
		g = ref_genomes[gname]
		
		if args.genome is not None and args.genome != gname:
			continue #continue if not specific genome
		
		data = []
		data.append({
			'Attribute': 'Genome Size',
			'Value': g.gsize,
			'Exists?': 'N/A'
		})
		data.append({
			'Attribute': 'Chrom Sizes',
			'Value': g.chrsize,
			'Exists?': str(os.path.exists(g.chrsize) and os.path.isfile(g.chrsize))
		})
		data.append({
			'Attribute': 'Whole Genome FASTA',
			'Value': g.wg_fasta,
			'Exists?': str(os.path.exists(g.wg_fasta) and os.path.isfile(g.wg_fasta))
		})
		for chrom_fasta in sorted(g.chr_fasta):
			data.append({
				'Attribute': '{} FASTA'.format(chrom_fasta),
				'Value': g.chr_fasta[chrom_fasta],
				'Exists?': str(os.path.exists(g.chr_fasta[chrom_fasta]) and os.path.isfile(g.chr_fasta[chrom_fasta]))
			})
		data.append({
			'Attribute': 'Gene Annotation File',
			'Value': g.genes_gtf,
			'Exists?': str(os.path.exists(g.genes_gtf) and os.path.isfile(g.genes_gtf))
		})
		for idx in sorted(g.indexes):
			data.append({
				'Attribute': idx,
				'Value': g.indexes[idx],
				'Exists?': str(len(glob.glob(g.indexes[idx]+'*')) > 0)
			})
		
		sys.stderr.write('Inspecting Genome "{}":\n'.format(gname))
		sys.stderr.write(tabulate.tabulate(data, headers='keys', tablefmt="grid"))
		sys.stderr.write('\n\n\n')
	







if __name__ == "__main__":
	main()










	
