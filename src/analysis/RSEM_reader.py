#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
# extract lines from RSEM exon data which fall within regions of interest as defined by a BED file
# * optionally convert first column of RSEM format into columns (chrom, start, end, gene, [strand]), with gene name from BED

# Usage: RSEM_reader.py [options] regions.bed
# -i infn:  Input data filename
# -o outfn:  Output data filename
# -d  Turn on debugging
# -c  Convert first column of RSEM data into columns [chr, start, end, gene]
# -p  Add strand information to leading columns.  Implies -c
# -D  Exclude data columns, write only leading columns.
# -H  Do not write header
# -s N:  Remove N leading characters in chrom name

import sys

# Sample BED:
# 14  68649185    68792100
class BED_filter:
    """Reads a bed file consisting of one or more lines (segments) with chromosome, start, end, and optional name.
    Tests whether a given position and chromosome is in this region of interest.
    """
    def __init__(self, bedfn):
        f = open(bedfn, "r")
        self.segments = []   # Segments is a list of (chr, begin, end, name) tuples, one for each line in BED file
        for l in f:
            tok = l.rstrip().split("\t")
            c, b, e = tok[0], int(tok[1]), int(tok[2])
            s = tok[3] if (len(tok) > 3) else ''
            self.segments.append( (c,b,e,s) )
        f.close()

    def inROI(self, seg, c, pos):
        # is given chromosome and position in region of interest?
        seg_chr, seg_begin, seg_end, seg_name = seg
        if (not c  == seg_chr): return False
        if (pos < seg_begin): return False
        if (pos > seg_end): return False
        return True

    # A line is printable if it is in region of interest (i.e., in proper chr and between start and end positions)
    def printable(self, c, p):
        return False if self.get_segment(c,p) == None else True

    # Get first BED file segment within which this chrom/pos falls.  
    # Segment is tuple (chrom, start, stop, name)
    # Return None if no matching segment found
    def get_segment(self, c, p):
        # c: chromosome of depth line
        # p: genome position of depth line
        for seg in self.segments:
            if self.inROI(seg, c, p): return seg
        return None

def filter_exon(f, o, bf, options):
    for i, line in enumerate(f, start=1):
    # Sample RSEM line:
    # chr14:68290259-68290344:+
        col1, data = line.rstrip().split("\t", 1)

        if options.no_data:
            data = ""

        # Skip header row starting with "exon"
        if (col1 == "exon"): continue

        # Header rows start with "Hybridization REF"
        # Print as necessary and continue to next line
        if (col1 == "Hybridization REF"): 
            if options.strand:
                col1 = '\t'.join( ('chrom', 'start', 'end', 'gene', 'strand') )
            elif options.convert_columns:
                col1 = '\t'.join( ('chrom', 'start', 'end', 'gene') )
            if not options.no_header:
                o.write('\t'.join( (col1, data)) + "\n")
            continue

        # Process the first column 
        tok1 = col1.split(":")
        if len(tok1) != 3:
            raise Exception("Weird column 1 format: %s" % t)

        (chrom, pos, strand) = tok1

        pos_start, pos_end = map(int, pos.split("-"))
        mid = pos_start + int((pos_end - pos_start) / 2)
        if options.strip_chr:
            chrom = chrom[int(options.strip_chr):] 

        if options.debug: print "Testing chr %s, %d - %d" % (chrom, pos_start, pos_end)

        matched_segment = bf.get_segment(chrom, mid)
        if matched_segment is not None:
            matched_chrom, matched_start, matched_end, matched_gene = matched_segment
            if options.strand:
                col1 = '\t'.join( map(str, (chrom, pos_start, pos_end, matched_gene, strand)) )
            elif options.convert_columns:
                col1 = '\t'.join( map(str, (chrom, pos_start, pos_end, matched_gene)) )

            o.write('\t'.join( (col1, data)) + "\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] BED.fn ...
        Extract lines of RSEM exon file where exon overlaps with
        domains specified in BED file. Optionally expand first column 
        of RSEM file into RPKM format, with columns [chr, start, end, gene],
        where gene name from BED file. 

        By default reads from stdin, writes to stdout (or -i, -o)
       
        Example: tar -zxOf origdata/HNSC_exon.tar.gz | python RSEM_reader.py BED.fn > dat/BA-4077_RSEM.dat
        
        """

    parser = OptionParser(usage_text)
    parser.add_option("-i", dest="infn", default="stdin", help="Input data filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output data filename")

    parser.add_option("-d", dest="debug", default=False, action="store_true", help="Turn on debugging")
    parser.add_option("-c", dest="convert_columns", default=False, action="store_true", help="Convert first column of RSEM data into columns [chr, start, end, gene]")
    parser.add_option("-p", dest="strand", default=False, action="store_true", help="Add strand information to leading columns.  Implies -c") 
    parser.add_option("-D", dest="no_data", default=False, action="store_true", help="Exclude data columns, write only leading columns.")  
    parser.add_option("-H", dest="no_header", default=False, action="store_true", help="Do not write header")  
    parser.add_option("-s", dest="strip_chr", help="Remove N leading characters in chrom name")  

    (options, params) = parser.parse_args()

    if (len(params) != 1):
        parser.error("Pass one arguments.")
    bedfn = params[0]

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    bf = BED_filter(bedfn)


    filter_exon(f,o,bf,options)



if __name__ == '__main__':
    main()
