#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu

import sys

class ChromTranslator:
    """ Normalize chrom names according to map file """

    def __init__(self, chrom_names_fn, key_col_second):
        self.chrom_map = {}
        for line in open(chrom_names_fn, 'r'):
            # 1   chr1
            t = line.rstrip().split("\t")
            (k,v) = (t[0],t[1]) if not key_col_second else (t[1], t[0])
            self.chrom_map[k] = v

    def rename_chrom(self, chrom):
        """ if chrom is a key, return corresponding value, otherwise return key """
        if chrom in self.chrom_map.keys():
            return self.chrom_map[chrom]
        else:
            return chrom

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] chrom_names.dat
        Translate chromosome names in BED file between two standards using chrom_names as database:
            names from first column of input data (-i) are replaced by translated chrom name;
            any instance of name from first column of database replaced by value from second column.

        GRCh38 Ensembl to Gencode conversion can be obtained here:
            https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2gencode.txt
        """

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-B", dest="key_col_second", action="store_true", help="Key name in second column of chrom_names.dat (first otherwise)")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename, BED format")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename, BED format")

    (options, params) = parser.parse_args()

    if (len(params) != 1):
        parser.error("Pass 1 argument.")

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    tr = ChromTranslator(params[0], options.key_col_second)

    # Parse BED file
    for line in f:
        # Skip headers
        if line.startswith("#"): next

        t = line.rstrip().split("\t")

        o.write( '\t'.join( ([tr.rename_chrom(t[0])] + t[1:] + ["\n"]) ))

    f.close()
    o.close()

if __name__ == '__main__':
    main()

