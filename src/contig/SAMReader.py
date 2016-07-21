#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu

# TODO: make the output of this be 1-base rather than 0-base, so that this output is consistent with
# what is seen in SAM file.  This will make SBP be 1-base too, making analysis easier.

# Usage:
#   python BreakpointReader.py [-o out.dat] data.sam
# Version 1.6 5/12/15  - drops requirement of only one breakpoint per contig
# Version 1.0 5/3/15  - first working version

# Extract breakpoint positions from SAM files with split reads
# All coordinates are printed in the native coordinates of pysam (0 index)

import sys
import pysam
# http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment    

# CIGAR codes
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8

class pSBPwriter:
    """Responsible for output formatting"""
    # The following tab-delimited fields are written to the file o,
    # QNAME, RNAME, SEQ_LEN, START_SPOS, END_SPOS, START_RPOS, END_RPOS, IS_REVERSE, IS_PRIMARY
    # where the first two fields are as in SAM file. 
    # SEQ_LEN is the query length. It is most meaningful for IS_PRIMARY=TRUE
    # START_SPOS, END_SPOS are the start, end sequence positions as given by CIGAR string for this alignment
    # and  START_RPOS, END_RPOS are corresponding reference positions 
    # IS_REVERSE indicates the status of the reverse complement flag (0x10)
    # IS_PRIMARY is NOT(IS_SUPPLEMENTARY), given by flag 0x800

    header = ("query_name", "ref_name", "seq_len", "start_spos", "end_spos", "start_rpos", "end_rpos", "is_reverse", "is_primary")

    def __init__(self, o):
        """Instantiate writer by passing output object"""
        self.o = o

    def writeLine( self, qname, rname, seq_len, start_spos, end_spos, start_rpos, end_rpos, is_reverse, is_primary ):
        # pass all values necessary to write a line and write to output device
        data = ( qname, rname, seq_len, start_spos, end_spos, start_rpos, end_rpos, is_reverse, is_primary )
        self.o.write('\t'.join( map(str, data))+"\n")

    def writeHeader(self):
        self.writeLine(*self.header)  # http://stackoverflow.com/questions/1993727/expanding-tuples-into-arguments

def isClipped(sam):
    """Returns true if any segment in CIGAR is hard or soft clipped"""
#    for (op,l) in sam.cigartuples:
#        print op, l
    areClipped = [op in (BAM_CSOFT_CLIP, BAM_CHARD_CLIP) for (op,l) in sam.cigartuples]
    return any(areClipped)

def getLeadingHardClipPadding(sam):
    """Return length of leading hard clipped segment, 0 if no hard clipping"""
    # hard clipping is not reflected in get_reference_positions(), so get leading padding
    # to give correct position of alignment in contig
    return sam.cigartuples[0][1] if sam.cigartuples[0][0] == BAM_CHARD_CLIP else 0

def parseSAM(samFn, w, mapq) :
    # for each line of the sam file corresponding to an alignment of a contig, we evaluate the start and end reference as well
    # as sequence positions.  Alignments with a mapping quality less than that of 'mapq' filter are discarded.
    # Alignments with no clipping in CIGAR string (implying no breakpoints) are also discarded.
    # Note that sequence position is in "local" coordinates with respect to alignment of this segment; segments which have
    # a direction (as given by is_reverse) which is opposite to the primary segment (!is_supplementary) will have start, end sequence
    # positions swapped.
    # Information from one line of pSBP is not enough to yield breakpoint position; this can only be determined downstream
    # when a pair of segments are combined.
    # http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment
    samfile = pysam.Samfile( samFn, "r" )

    #process all reads in region
    for sam in samfile.fetch():
        if sam.is_unmapped: 
            continue
        if not isClipped(sam):
            continue
        if sam.mapping_quality < mapq:
            continue

        refPosMapped = sam.get_reference_positions()
        refPosAll = sam.get_reference_positions(full_length=True)   # has None for soft clipped regions
        
        start_rpos = refPosMapped[0]  # first mapped reference position
        end_rpos = refPosMapped[-1]  # last mapped reference position
        # get first, last mapped sequence positions by an inverse lookup, and add leading hard clip
        # http://stackoverflow.com/questions/176918/finding-the-index-of-an-item-given-a-list-containing-it-in-python
        start_spos = refPosAll.index(start_rpos) + getLeadingHardClipPadding(sam)
        end_spos = refPosAll.index(end_rpos) + getLeadingHardClipPadding(sam)

        seq_len = sam.query_length   # most meaningful for primary alignments 

        rname = samfile.getrname(sam.rname) 
        w.writeLine( sam.query_name, rname, seq_len, start_spos, end_spos, start_rpos, end_rpos, sam.is_reverse, not sam.is_supplementary)

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] data.sam ...
        Extract breakpoint positions from a SAM file of a chimeric contig aligned to multiple references
        """

    parser = OptionParser(usage_text, version="$Revision: 1.0 $")
    parser.add_option("-v", dest="verbose", help="Print extra output")
    parser.add_option("-o", dest="out_fn", default="stdout", help="Write output in TSV format [stdout]")
    parser.add_option("-Q", dest="MAPQfilter", default="60", help="Discard alignments with MAPQ score < given value")

    (options, params) = parser.parse_args()

    if options.verbose:
        print options
        print params
    mapq = int(options.MAPQfilter)

    if (len(params) != 1):
        parser.error("Pass 1 argument.")
    dataFn = params[0]

    if ( options.out_fn == "stdout" ):
        o = sys.stdout
    else:
        o = open(options.out_fn, 'w')
        print "Writing to ", options.out_fn
    w = pSBPwriter(o)
    w.writeHeader()
    parseSAM(dataFn, w, mapq)

if __name__ == '__main__':
    main()
