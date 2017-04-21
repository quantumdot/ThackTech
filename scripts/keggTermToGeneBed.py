#!/usr/bin/env python

import argparse
import MySQLdb
import sys

def main():
    

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', required=True, help="Reference genome build to work with, in UCSC terms (i.e. mm9).")
    parser.add_argument('-t', '--term', required=True, help="KEGG pathway term to search for. By default looks for an exact match.")
    parser.add_argument('outfile', help="File to write results in BED format")
    
    db_conn_group = parser.add_argument_group("SQL Connection Options")
    db_conn_group.add_argument('--host', default='genome-mysql.cse.ucsc.edu', help="hostname of the SQL server")
    db_conn_group.add_argument('--user', default='genome', help="username of the SQL server")
    db_conn_group.add_argument('--port', default=3306, type=int, help="username of the SQL server")
    args = parser.parse_args()
    
    
    connection = MySQLdb.connect(host=args.host, user=args.user, port=args.port, db=args.genome)
    
    try:
        sub_query = "SELECT keggPathway.kgID " \
                  + "FROM keggPathway INNER JOIN keggMapDesc ON keggMapDesc.mapID = keggPathway.mapID " \
                  + "WHERE keggMapDesc.description = '{query}'".format(query=connection.escape_string(args.term))
                  
        sql = "SELECT chrom, txStart, txEnd, CONCAT(kgXref.mRNA, '|', kgXref.refseq, '|', kgXref.geneSymbol) as name, '.' as score, strand " \
            + "FROM knownGene  INNER JOIN kgXref ON knownGene.name = kgXref.kgID " \
            + "WHERE kgXref.kgID IN ({subquery})".format(subquery=sub_query)
        
        cursor = connection.cursor()
        cursor.execute(sql)
        
        if cursor.rowcount > 0:
            sys.stderr.write("Found {num_results} results\nWriting to {path}....\n".format(num_results=cursor.rowcount, path=args.outfile))
            with open(args.outfile, 'w') as outfile:
                for data in cursor.fetchall():
                    outfile.write("\t".join([str(d) for d in data]))
                    outfile.write("\n")
        else:
            sys.stderr.write("WARNING: No results were returned for query term '{query}'!!!\n".format(query=args.term))
    finally:
        connection.close()
    
    
    
if __name__ == "__main__":
    main()