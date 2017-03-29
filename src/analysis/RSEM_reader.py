#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
# extract lines from RSEM exon data which fall within regions of interest as defined by a BED file
# Optionally convert first column of RSEM format into columns [chrom, start, end, gene], with gene name from BED

# v2.1: optionally convert first column of RSEM format into columns [chrom, start, end, gene], with gene name from BED
# v2.0: retains exon data based on domains defined in BED file.

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
        t, data = line.rstrip().split("\t", 1)

        # Skip header row starting with "exon"
        if (t == "exon"): continue

        # retain and expand if necessary header rows starting with "Hybridization REF"
        if (t == "Hybridization REF"): 
            if options.convert_columns:
                t = '\t'.join( ('chr', 'start', 'end', 'gene') )
            o.write('\t'.join( (t, data)) + "\n")
            continue

        # Process the first column entry
        tok1 = t.split(":")
        if len(tok1) != 3:
            raise Exception("Weird column 1 format: %s" % t)
        (c, r, p) = tok1
        RSEMstart, RSEMend = map(int, r.split("-"))
        pos = RSEMstart + int((RSEMend - RSEMstart) / 2)
        RSEMchr = c[3:] if options.strip_chr else c
        if options.debug: print "Testing chr %s, %d - %d" % (RSEMchr, RSEMstart, RSEMend)

        seg = bf.get_segment(RSEMchr, pos)
        if seg is not None:
            if options.convert_columns:
                t = '\t'.join( map(str, (RSEMchr, RSEMstart, RSEMend, seg[3])) )

            o.write('\t'.join( (t, data)) + "\n")

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

    #  TODO: Add -i, -o to input/output filenames

    parser = OptionParser(usage_text, version="$Revision: 2.0 $")
    parser.add_option("-d", dest="debug", default=False, action="store_true", help="Turn on debugging")
    parser.add_option("-c", dest="convert_columns", default=False, action="store_true", help="Convert first column of RSEM data into columns [chr, start, end, gene]")
    parser.add_option("-s", dest="strip_chr", default=False, action="store_true", help="Remove leading 'chr' in chrom name")
    parser.add_option("-i", dest="infn", default="stdin", help="Input data filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output data filename")


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
