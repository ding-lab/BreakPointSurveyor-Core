#!/usr/bin/python
# Matthew A. Wyczalkowski, m.wyczalkowski@wustl.edu

# Create a Breakpoint Surveyor PlotList file from Breakpoint Coordinate (BPC) or Breakpoint Region (BPR) data
# Each line in a PlotList corresponds to one plot to be generated.  
# Breakpoints correspond to breakpoint coordinates A, B, with A and B typically on different chromosome 
#  * A and B chrom can be the same, and either can be a virus 
#  * by default chrom A < chrom B in a string comparison sense; this is reversed with -l flag
# Plots are generated around some events of interest (e.g. SV regions), and plot ranges are given by event regions +/- "context" distance
# PlotList is TSV format with the following columns,
#  * barcode
#  * event.name (unique)
#  * chrom.a, (first chromosome of coordinate pair)
#  * event.a.start, event.a.end (indicates region of e.g. SV event)
#  * range.a.start, range.a.end (indicates region to plot; calculated as event.start - context, event.end + context, respectively)
#  * chrom.b, (second chromosome of coordiante pair)
#  * event.b.start, event.b.end, range.b.start, range.b.end 
# attribute columns ignored 

# BPC, BPR data formats: $BPS_CORE/doc/BPC_BPR_FileFormat.txt

import itertools, string, sys

class FAI:
# reads, stores, and processes chromosome lengths based on FAI file.
# sample line: 1   249250621   52  60  61
# where first col is chrom name, second is chrom length.  Remaining ignored.

    def __init__(self, fai_fn):
        """Reads in FAI file.  If fai_fn is None, subsequent cropping returns values passed"""

        if fai_fn is not None:
            self.chrom_lengths = {}
            f = open(fai_fn, 'r')
            for l in f: 
                (chrom, chrom_len, a, b, c) = l.split("\t")
                self.chrom_lengths[chrom] = int(chrom_len)
            f.close()
        else:
            self.chrom_lengths = None

    def crop(self, chrom, start, end):
        """Return start, end positions of given chrom which correspond to valid genomic coordinates.
        
        Require that start >= 1 and end <= chromosome_length """
        if self.chrom_lengths is None:
            return start, end
        else:
            return max(1, start), min(end, self.chrom_lengths[chrom])


def parse_BP(f, o, fai, isBPC, barcode, context, options):
    # suffix is series of strings "AA", "AB", "AC", etc. which ensure event names are unique
    # a generator and http://stackoverflow.com/questions/2156892/python-how-can-i-increment-a-char
    suffix=(''.join(i) for i in itertools.product(string.ascii_uppercase, repeat=int(options.suffix_length)))  

    if options.write_header:
        o.write('\t'.join( ("barcode", "name", "chrom.A", "event.A.start", "event.A.end", "range.A.start", "range.A.end",
                                                   "chrom.B", "event.B.start", "event.B.end", "range.B.start", "range.B.end") ))
        o.write('\n')
    for line in f:
        if line[0] == "#": continue
        t = line.rstrip().split("\t")
        if t[0] == "chromA": continue  # header

#     BPC: chromA, posA, chromB, posB, [attrib]
#     BPR: chromA, posA.start, posA.end, chromB, posB.start, posB.end, [attrib]
        if isBPC:
            if len(t) == 4:
                chromA, posA_start, chromB, posB_start = t
            elif len(t) == 5:
                chromA, posA_start, chromB, posB_start, attrib = t
            else:
                raise Exception("Unknown format of BPC input file")
            posA_start, posB_start = map(int, (posA_start, posB_start))
            posA_end, posB_end = posA_start + 1, posB_start+1
        else:
            if len(t) == 6:
                chromA, posA_start, posA_end, chromB, posB_start, posB_end = t
            elif len(t) == 7:
                chromA, posA_start, posA_end, chromB, posB_start, posB_end, attrib = t
            else:
                raise Exception("Unknown format of BPR input file")
            posA_start, posA_end, posB_start, posB_end = map(int, (posA_start, posA_end, posB_start, posB_end))

#  * event.name (unique)
#  * barcode
#  * chrom.a, (first chromosome of coordinate pair)
#  * event.a.start, event.a.end (indicates region of e.g. SV event)
#  * range.a.start, range.a.end (indicates region to plot; calculated as event.start - context, event.end + context, respectively)
#  * chrom.b, (second chromosome of coordiante pair)
#  * event.b.start, event.b.end, range.b.start, range.b.end 
        try:
            # chrom_keep defines which chrom names listed in event name.  May exclude names if they are long (e.g., virus) or otherwise unnecessary
            chrom_keep = []
            if "A" in options.chrom_keep: chrom_keep.append(chromA)
            if "B" in options.chrom_keep: chrom_keep.append(chromB)
            # chrom_prefix can be "chr" if chrom name is e.g. "1".  This makes event name format more clear.
            chrom_names = [options.chrom_prefix + c for c in chrom_keep]
            event_name = ".".join( (barcode, suffix.next(), "_".join( chrom_names )) )
        except StopIteration: # Catch StopIteration, require suffix_length to be increased
            sys.stderr.write("Event name looped.  Increase suffix_length (-s)\n")
            sys.exit()

        rA_start, rA_end = fai.crop(chromA, posA_start - context, posA_end + context)
        rB_start, rB_end = fai.crop(chromB, posB_start - context, posB_end + context)

        # Default order of columns is given by BPR/BPC convention: chromA < chromB using alphabetical comparison.
        # flip_ab (-l) reverses this order
        # XOR operator: True ^ True = False, True ^ False = True, False ^ True = True, False ^ False = False
        if (str(chromA) == str(chromB)):
            a_then_b = (posA_start < posB_start) ^ options.flip_ab
        else:
            a_then_b = (str(chromA) < str(chromB))  ^ options.flip_ab

        if a_then_b:
            cols = (barcode, event_name, chromA, posA_start, posA_end, rA_start, rA_end, chromB, posB_start, posB_end, rB_start, rB_end)
        else:
            cols = (barcode, event_name, chromB, posB_start, posB_end, rB_start, rB_end, chromA, posA_start, posA_end, rA_start, rA_end)
        o.write('\t'.join( map(str, cols) ) )
        o.write('\n')

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] ...
        Generate a PlotList file from a BPS or BPR file

        Unique names are created for each line of PlotList file as,
            BARCODE.SUFFIX_CPA_CPB,
        where SUFFIX is AA, AB, ...; CPA and CPB are composed of chrom_prefix+chrom_name
        """

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-c", dest="context", default="50000", help="Padding around breakpoint regions")
    parser.add_option("-C", dest="BPC", action="store_true", help="Input file is BPC format (default BPR)")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename")
    parser.add_option("-r", dest="reference_fai", default=None, help="Reference FAI file listing chrom lengths")
    parser.add_option("-n", dest="barcode", default="", help="Sample barcode or other identifier")
    parser.add_option("-s", dest="suffix_length", default="2", help="Length of letter code (2 is aa, ab, ac, ...)")
    parser.add_option("-p", dest="chrom_prefix", default="", help="String before chrom name in event name, for legibility")
    parser.add_option("-H", dest="write_header", action="store_true", help="Write column headers in output")
    parser.add_option("-N", dest="chrom_keep", default="AB", help="One of AB, A, B, X; which chrom names to include in name, X for none")
    parser.add_option("-l", dest="flip_ab", action="store_true", default=False, help="Swap chrom A and B ")

    (options, params) = parser.parse_args()

    if options.infn == "stdin":
        f = sys.stdin
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    # this works fine even if reference_fai not passed
    fai = FAI(options.reference_fai)

    parse_BP(f, o, fai, options.BPC, options.barcode, int(options.context), options)

    f.close()
    o.close()

if __name__ == '__main__':
    main()

