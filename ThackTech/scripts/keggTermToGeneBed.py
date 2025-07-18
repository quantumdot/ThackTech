#!/usr/bin/env python
import os
import sys
import urllib
import urllib2
import argparse
import MySQLdb
import MySQLdb.cursors

from ThackTech import conf
config = conf.get_config('ontology_to_genes')



def main():
    

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', required=True, choices=['go', 'kegg'], help="annotation database to search")
    parser.add_argument('-g', '--genome', required=True, help="Reference genome build to work with, in UCSC terms (i.e. mm9).")
    parser.add_argument('-t', '--term', required=True, help="KEGG pathway term to search for. By default looks for an exact match.")
    parser.add_argument('--type', choices=['name', 'acc'], default='acc', help='Is the query a term name or accession number')
    parser.add_argument('--gomode', choices=['direct', 'indirect'], default='direct', help="For GO queries, return directly annotated genes or include transitive indirect annotated genes.")
    parser.add_argument('outfile', help="File to write results in BED format")
    
    #db_conn_group = parser.add_argument_group("SQL Connection Options")
    #db_conn_group.add_argument('--host', default='genome-mysql.cse.ucsc.edu', help="hostname of the SQL server")
    #db_conn_group.add_argument('--user', default='genome', help="username of the SQL server")
    #db_conn_group.add_argument('--port', default=3306, type=int, help="username of the SQL server")
    args = parser.parse_args()
    
    if args.database == 'go':
        results = genes_for_go_term(args.term, args)
    elif args.database == 'kegg':
        results = genes_for_kegg_term(args.term, args)
        
    results = filter_genes_canonical(results)
    
    if results is not None and len(results) > 0:
        sys.stderr.write("Found {num_results} results\nWriting to {path}....\n".format(num_results=len(results), path=args.outfile))
        with open(args.outfile, 'w') as outfile:
            for data in results:
                outfile.write("\t".join([str(d) for d in data]))
                outfile.write("\n")
    else:
        sys.stderr.write("WARNING: No results were returned for query term '{query}'!!!\n".format(query=args.term))
    
#end main()

__connections = {}
def get_connection(name):
    if name not in __connections:
        options = {
            'host': config.get(name, 'host'),
            'user': config.get(name, 'user'),
            'passwd': config.get(name, 'pass'),
            'port': config.getint(name, 'port'),
            #'cursorclass': MySQLdb.cursors.DictCursor
        }
        c = MySQLdb.connect(**options)
        __connections[name] = c
        
    return __connections[name]
#end get_connection_opts()

def fetch_results(connection_name, database, query, data=None):
    db = get_connection(connection_name)
    db.select_db(database)
    
    cursor = db.cursor()
    cursor.execute(query, data)
    
    if cursor.rowcount > 0:
        return cursor.fetchall()
    else:
        return []
#end fetch_results()

def ucsc_genome_to_ncbi_taxid(genome):
    db = get_connection('ucsc_connection')
    db.select_db('hgcentral')
    sql = "SELECT DISTINCT taxId FROM dbDb WHERE name = '{query}'".format(query=db.escape_string(genome))
    cursor = db.cursor()
    cursor.execute(sql)
    if cursor.rowcount > 0:
        return cursor.fetchone()[0]
    else:
        return None
#end ucsc_genome_to_ncbi_taxid()

def genes_for_kegg_term(term, options):
    sub_query = "SELECT keggPathway.kgID " \
              + "FROM keggPathway INNER JOIN keggMapDesc ON keggMapDesc.mapID = keggPathway.mapID " \
              + "WHERE keggMapDesc.description = '%s'"
              
    sql = "SELECT chrom, txStart, txEnd, CONCAT(kgXref.mRNA, '|', kgXref.refseq, '|', kgXref.geneSymbol) as name, '.' as score, strand " \
        + "FROM knownGene INNER JOIN kgXref ON knownGene.name = kgXref.kgID " \
        + "WHERE kgXref.kgID IN ("+sub_query+")"
    
    return fetch_results('ucsc_connection', options.genome, sql, [term])
#end genes_for_kegg_term()

def genes_for_refseq_ids(refseq_ids, options):   
    format_strings = ','.join(['%s'] * len(refseq_ids)) 
    #sql = "SELECT chrom, txStart, txEnd, CONCAT(kgXref.mRNA, '|', kgXref.refseq, '|', kgXref.geneSymbol) as name, '.' as score, strand " \
    #    + "FROM knownGene INNER JOIN kgXref ON knownGene.name = kgXref.kgID " \
    #    + "WHERE kgXref.refseq IN (%s)" % (format_strings,)
    sql = "SELECT chrom, txStart, txEnd, CONCAT(name, '|', geneName) as name, '.' as score, strand " \
        + "FROM refFlat " \
        + "WHERE name IN (%s)" % (format_strings,)
    return fetch_results('ucsc_connection', options.genome, sql, tuple(refseq_ids))
#end genes_for_refseq_ids()

def filter_genes_canonical(bedlike_entries):
    genes_by_id = {}
    max_gene_lens = {}
    
    for b in bedlike_entries:
        name = b[3].split('|')[1]
        if name not in genes_by_id:
            genes_by_id[name] = []
            max_gene_lens[name] = 0
        max_gene_lens[name] = max(max_gene_lens[name], b[2] - b[1])
        genes_by_id[name].append(b)
        
    results = []
    for group in genes_by_id:
        for locus in genes_by_id[group]:
            if locus[2] - locus[1] == max_gene_lens[group]:
                results.append(locus)
                break
    sort_genes(results)
    return results
#end filter_genes_canonical()

def sort_genes(bedlike_entries):
    sorted(bedlike_entries, key=lambda e: (e[0], e[1]))
#end sort_genes
    
__id_mapping = {
    'BioCyc': 'BIOCYC_ID',
    'EMBL': 'EMBL_ID',
    'Ensembl': 'ENSEMBL_ID',
    'GeneDB': 'GENEDB_ID',
    'MGI': 'MGI_ID',
    'UniProtKB': 'ID',
    'PDB': 'PDB_ID'
}
def convert_ids_to_refseq(ids_by_source):
    """Expect a dict of lists with key = db source and list of ids
    """
    url = "http://www.uniprot.org/uploadlists/"
    results = []
    for source in ids_by_source:
        if source not in __id_mapping:
            sys.stderr.write("WARNING: Could not convert IDs from source {}".format(source))
        
        if source != 'UniProtKB':
            dest = 'ID'
        else:
            dest = 'REFSEQ_NT_ID'
        
        data = urllib.urlencode({
            'from': __id_mapping[source],
            'to': dest,
            'format':'list',
            'query': ' '.join(ids_by_source[source])
        })
        request = urllib2.Request(url, data)
        request.add_header('User-Agent', 'Python')
        response = urllib2.urlopen(request)
        
        if source != 'UniProtKB':
            results.extend(convert_ids_to_refseq({'UniProtKB': [line.strip() for line in response]}))
        else:
            results.extend([os.path.splitext(line.strip())[0] for line in response])
    return results
#end uniprot_to_refseq()



def get_go_query_include_transitive(search_field):
    """Gets a SQL query for a GO term including transitive indirect annotations
    
    Adapted from: http://wiki.geneontology.org/index.php/Example_LEAD_Queries#All_genes_in_Drosophila_annotated_to_.27nucleus.27_.28including_transitive_indirect_annotations.29
    """
    sql = "SELECT " \
        + "    dbxref.xref_dbname AS gp_dbname, " \
        + "    dbxref.xref_key AS gp_acc " \
        + "FROM term " \
        + "    INNER JOIN graph_path ON (term.id=graph_path.term1_id) " \
        + "    INNER JOIN association ON (graph_path.term2_id=association.term_id) " \
        + "    INNER JOIN gene_product ON (association.gene_product_id=gene_product.id) " \
        + "    INNER JOIN species ON (gene_product.species_id=species.id) " \
        + "    INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id) " \
        + "WHERE " \
        + "    term."+search_field+" = %s " \
        + "AND species.ncbi_taxa_id = %s"
    return sql
#end get_go_query_include_transitive()

def get_go_query_direct(search_field):
    """Gets a SQL query for a GO term including transitive indirect annotations
    
    Adapted from: http://wiki.geneontology.org/index.php/Example_LEAD_Queries#All_genes_directly_annotated_to_.27nucleus.27_.28excluding_child_terms.29
    """
    sql = "SELECT " \
        + "    dbxref.xref_dbname AS gp_dbname, " \
        + "    dbxref.xref_key AS gp_acc " \
        + "FROM term " \
        + "    INNER JOIN association ON (term.id=association.term_id) " \
        + "    INNER JOIN gene_product ON (association.gene_product_id=gene_product.id) " \
        + "    INNER JOIN species ON (gene_product.species_id=species.id) " \
        + "    INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id) " \
        + "    INNER JOIN db ON (association.source_db_id=db.id) " \
        + "WHERE " \
        + "    term."+search_field+" = %s " \
        + "AND species.ncbi_taxa_id = %s"
    return sql
#end get_go_query_direct()

def genes_for_go_term(term, options):
    taxid = ucsc_genome_to_ncbi_taxid(options.genome)
    
    if options.gomode == 'indirect':
        sql = get_go_query_include_transitive(options.type)
    else:
        sql = get_go_query_direct(options.type)
    
    #sys.stderr.write(sql+"\n")
    go_hits = set(fetch_results('go_connection', 'go_latest', sql, (term, taxid)))
    if len(go_hits) == 0:
        sys.stderr.write("Did not find any GO matches for term {}!\nPlease check that it exists!\n".format(term))
        return []
    ids_by_source = {}
    for h in go_hits:
        if h[0] not in ids_by_source:
            ids_by_source[h[0]] = []
        ids_by_source[h[0]].append(h[1])
        
    sys.stderr.write("{} go hits\n".format(len(go_hits)))
    ref_ids = set(convert_ids_to_refseq(ids_by_source))
    sys.stderr.write("{} refseq ids\n".format(len(ref_ids)))
    #print ref_ids
    return genes_for_refseq_ids(ref_ids, options)
    
if __name__ == "__main__":
    main()