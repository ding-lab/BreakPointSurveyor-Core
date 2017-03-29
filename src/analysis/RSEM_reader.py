#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
# extract lines from RSEM exon data which fall within regions of interest as defined by a BED file

# v2.0: retains exon data based on domains defined in BED file.

# BED filter adapted from, /gscuser/mwyczalk/projects/Virus/Virus_2013.9a/analysis/UnifiedVirus2/Q2_Breakpoints/src/depthFilter.py
import sys

# Sample BED:
# 14  68649185    68792100
class BED_filter:
    """Reads a bed file consisting of one or more lines (segments) with chromosome, start and end.
    Tests whether a given position and chromosome is in this region of interest.
    """
    def __init__(self, bedfn):
        f = open(bedfn, "r")
        self.segments = []   # Segments is a list of (chr, begin, end) tuples, one for each line in BED file
        for l in f:
            tok = l.rstrip().split("\t")
            c, b, e = tok[0], int(tok[1]), int(tok[2])
            self.segments.append( (c,b,e) )
        f.close()

    def inROI(self, seg, c, pos):
        # is given chromosome and position in region of interest?
        seg_chr, seg_begin, seg_end = seg
        if (not c  == seg_chr): return False
        if (pos < seg_begin): return False
        if (pos > seg_end): return False
        return True

    # A line is printable if it is in region of interest (i.e., in proper chr and between start and end positions)
    def printable(self, c, p):
        # c: chromosome of depth line
        # p: genome position of depth line
        for seg in self.segments:
            if self.inROI(seg, c, p): return True
        return False

def filter_exon(f, o, bf, options):
    for i, line in enumerate(f, start=1):
    # Sample RSEM line:
    # chr14:68290259-68290344:+
        t = line.rstrip().split("\t")[0]
        # Skip headers
        if (t == "Hybridization REF") or (t == "exon"): continue
        tok1 = t.split(":")
        if len(tok1) != 3:
            raise Exception("Weird header line %d: %s" % (i, t))
        (c, r, p) = tok1
        RSEMstart, RSEMend = map(int, r.split("-"))
        pos = RSEMstart + int((RSEMend - RSEMstart) / 2)
        RSEMchr = c[3:]
        if options.debug: print "Testing chr %s, %d - %d" % (RSEMchr, RSEMstart, RSEMend)
        if bf.printable(RSEMchr, pos): o.write(line)

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] BED.fn ...
        Extract lines of RSEM exon file where exon midpoint overlaps with
        domains specified in BED file.  
        Reads from stdin, writes to stdout.
       
        Example: tar -zxOf origdata/HNSC_exon.tar.gz | python RSEM_reader.py BED.fn > dat/BA-4077_RSEM.dat
        
        """

    parser = OptionParser(usage_text, version="$Revision: 2.0 $")
    parser.add_option("-d", dest="debug", default=False, action="store_true", help="Padding beyond break points where exon data is retained")

    (options, params) = parser.parse_args()

    if (len(params) != 1):
        parser.error("Pass one arguments.")
    bedfn = params[0]

    f = sys.stdin
    o = sys.stdout
    bf = BED_filter(bedfn)

    filter_exon(f,o,bf,options)



if __name__ == '__main__':
    main()
