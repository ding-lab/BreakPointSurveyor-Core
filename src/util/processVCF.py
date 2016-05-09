#!/usr/bin/python

# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute

import vcf
import sys

def write_size(o, f):
    """Write variant start, end, size, and type"""
    o.write( '\t'.join( ("Start", "End", "Size", "Type") ) + "\n")
    vcf_reader = vcf.Reader(f)
    for record in vcf_reader:
        o.write('\t'.join((str(record.start), str(record.end), str(record.end - record.start), record.var_subtype)) + "\n")

def write_bed(o, f):
    """Write variant BED file with variant chrom, start, end, and type"""
    o.write( "# " + '\t'.join( ("Chrom", "Start", "End", "Type") ) + "\n")
    vcf_reader = vcf.Reader(f)
    for record in vcf_reader:
        if isinstance(record.var_subtype, str):
            subtype = record.var_subtype  # sometimes is string
        else:
            subtype = record.var_subtype[0]  # sometimes is tuple
            # catch faulty assumptions.    
            if len(record.var_subtype) > 1: 
                print "Warning: multiple subtypes: ", record.var_subtype

        linedata = record.CHROM, record.start, record.end, subtype
        o.write('\t'.join( map(str, linedata) ) + "\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] mode ...
        mode is one of "bed", "size"
        """

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename")

    (options, params) = parser.parse_args()

    if (len(params) != 1):
        parser.error("Pass 1 argument.")
    mode = params[0]
    if mode not in ['bed', 'size']:
        parser.error("mode must be one of ('bed', 'size')")

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    o.write( "# " + options.infn + "\n" )
    if mode == "bed":
        write_bed(o, f)
    elif mode == "size":
        write_size(o, f)
            

    f.close()
    o.close()

if __name__ == '__main__':
    main()

