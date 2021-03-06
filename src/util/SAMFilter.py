#!/usr/bin/python
# Matthew A. Wyczalkowski, m.wyczalkowski@wustl.edu

# Based on contig/SAMReader.py and util/makeBreakpointRegions.py

### documentation ###

# Usage: SAMFilter.py [options] BPR.dat reads.sam ...
# Filter reads.sam to retain only those that fall within regions in BPR.dat
# 
# Options:
#   --version   show program's version number and exit
#   -h, --help  show this help message and exit


import sys
import pysam # http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment    

class BPR:
    """Read and store BPR data, and evaluate whether given point falls within it."""

    def __init__(self, fn, padding=0):
        self.BPR = []
        for line in open(fn, 'r'):
            if line.startswith("#"): continue
            self.BPR.append( self.parseBPR(line, padding) )

    def __repr__(self):
        return str(self.BPR)

    def parseBPR(self, line, padding=0):
        """Process one BPR line. Expand BPR by 'padding' amount, (start-padding, end+padding).
         Return (chromA, startA, endA, chromB, startB, endB, attr)
         If attribute not specified in file, attr=None """
        # Consider first seven columns, 
        t = line.rstrip().split("\t")
        attr = None if len(t) == 6 else t[6]
        padding = int(padding)
        return (t[0], int(t[1])-padding, int(t[2])+padding, t[3], int(t[4])-padding, int(t[5])+padding, attr)

    def contains(self, chromA, posA, chromB, posB):
        """Test whether chromA/posA and chromB/posB fall within any BPR region.  Assume A,B ordered in BPC format.  """
        # Note that if pos falls on a boundary it is considered within BPR
        # print "Testing " + str( (chromA, posA, chromB, posB))
        for (bpr_cA, bpr_sA, bpr_eA, bpr_cB, bpr_sB, bpr_eB, a) in self.BPR:
            # print "Against " + str( (bpr_cA, bpr_sA, bpr_eA, bpr_cB, bpr_sB, bpr_eB, a) )
            inA = (chromA == bpr_cA) and (posA >= bpr_sA) and (posA <= bpr_eA)
            inB = (chromB == bpr_cB) and (posB >= bpr_sB) and (posB <= bpr_eB)

            if inA and inB: return True

        return False


# See here for writing to/from stdin/out: http://pysam.readthedocs.io/en/latest/usage.html#using-streams


def filter_SAM(insam, outsam, bpr, options):

    #process all reads in region
    for read in insam.fetch():
        # We define position of read as first mapped position

        chromA = insam.get_reference_name(read.reference_id)
        posA = read.get_reference_positions()[0]
        chromB = read.next_reference_name
        posB = read.next_reference_start

        if bpr.contains(chromA, posA, chromB, posB):
            outsam.write(read)

# Return input and output SAM AlignmentFile objects 
def open_SAM(infn, outfn):

    if infn == "stdin": infn = "-"
    if outfn == "stdout": outfn = "-"

    insam = pysam.AlignmentFile(infn, "r")
    outsam = pysam.AlignmentFile(outfn, "w", template=insam)

    return (insam, outsam)

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] BPR.dat 
Filter reads to retain only those that fall within regions in BPR.dat
"""

# BPR-like because may have multiple attribute fields

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-d", dest="debug", action="store_true", help="Print cluster details as comments in output file")
    parser.add_option("-i", dest="infn", default="stdin", help="Input SAM file to parse.  Default reads from stdin")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output SAM file to write.  Default writes to stdout")
    parser.add_option("-p", dest="padding", default="0", help="Expand BPR regions by this amount")

    (options, params) = parser.parse_args()

    if (len(params) < 1):
        parser.error("Pass at least one argument.")
    bprfn = params[0]

    bpr = BPR(bprfn, options.padding)

    (insam, outsam) = open_SAM(options.infn, options.outfn)

    filter_SAM(insam, outsam, bpr, options)

    insam.close()
    outsam.close()


if __name__ == '__main__':
    main()
