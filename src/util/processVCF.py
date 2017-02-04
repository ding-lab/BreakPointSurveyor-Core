#!/usr/bin/python

# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute

# Read VCF file and write coordinates of features in various formats

import vcf # http://pyvcf.readthedocs.io/en/latest/API.html
import sys

def write_size(o, vcf_reader):
    """Write variant start, end, size, and type"""
    o.write( '\t'.join( ("Chrom", "Start", "End", "Size", "Type") ) + "\n")
    for record in vcf_reader:
        o.write('\t'.join((record.CHROM, str(record.start), str(record.end), str(record.end - record.start), record.var_subtype)) + "\n")

def write_bed(o, vcf_reader):
    """Write variant BED file with variant chrom, start, end, and type"""
    o.write( "# " + '\t'.join( ("Chrom", "Start", "End", "Type") ) + "\n")
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

def write_bpc(o, vcf_reader, attribute, no_header=False):
    """Write BreakpointSurveyor BPC file"""
#   From BreakpointSurveyor/doc/FileFormat.txt
#   BPC: chromA, posA, chromB, posB, [attribute]

# chr2    31823410    TRA00000050 N   <TRA>   .   PASS    IMPRECISE;SVTYPE=TRA;SVMETHOD=EMBL.DELLYv0.7.3;CHR2=chr1;END=33010827;INSLEN=0;PE=7;MAPQ=37;CT=5to3;CIPOS=-1009,1009;CIEND=-1009,1009   GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV   0/1:-17.3918,0,-3.81472:38:PASS:0:0:0:-1:1:7:0:0

    if not no_header:
        if attribute is None:
            o.write( "# " + '\t'.join( ("chromA", "posA", "chromB", "posB") ) + "\n")
        else:
            o.write( "# " + '\t'.join( ("chromA", "posA", "chromB", "posB", "attribute") ) + "\n")
        
    for record in vcf_reader:
        chromA, posA, chromB, posB = record.CHROM, record.start, record.INFO['CHR2'], record.sv_end
        if chromA < chromB:
            linedata = chromA, posA, chromB, posB 
        elif chromA > chromB:
            linedata = chromB, posB, chromA, posA 
        elif posA < posB:
            linedata = chromA, posA, chromB, posB 
        else:
            linedata = chromB, posB, chromA, posA 

        if attribute is not None:
            linedata += (attribute,)

        o.write('\t'.join( map(str, linedata) ) + "\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] mode ...
        mode is one of "bed", "size"
        """

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename")
    parser.add_option("-a", dest="attribute", default=None, help="Optional BPC attribute field")
    parser.add_option("-H", dest="no_header", action="store_true", help="Do not write header")

    (options, params) = parser.parse_args()

    if (len(params) != 1):
        parser.error("Pass 1 argument.")
    mode = params[0]
    if mode not in ['bed', 'size', 'bpc']:
        parser.error("mode must be one of ('bed', 'size', 'bpc')")

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    vcf_reader = vcf.Reader(f)
    if not options.no_header:
        o.write( "# " + options.infn + "\n" )
    if mode == "bed":
        write_bed(o, vcf_reader)
    elif mode == "size":
        write_size(o, vcf_reader)
    elif mode == "bpc":
        write_bpc(o, vcf_reader, options.attribute, no_header=options.no_header)
            

    f.close()
    o.close()

if __name__ == '__main__':
    main()

