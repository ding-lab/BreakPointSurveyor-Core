#!/usr/bin/python
# Matthew A. Wyczalkowski, m.wyczalkowski@wustl.edu

# Parse VAF as output by pindel
#
# TODO: make sure pindel, tigra parsers have common codebase.
#
# Example line:
#  1    5156542 .   .   .   60  PASS    IMPRECISE;SOMATIC;SVTYPE=TRA;CHR2=2;END=226658024;SVLEN=0

import sys

def filterVAF(f, o, options):
    for l in f:
        l = l.strip()
        if l[0] == "#": continue

        # chr1, pos1, lid, ref, alt, qual, lfilter, info, lformat, spikein, x = l.split("\t")  # this works for tigra-ext
        chr1, pos1, xa, xb, xc, xd, xe, info = l.split("\t")


        # split INFO string by ';', and create list of key/value tuples defined by =
        # create a dictionary out of those values which have length 2 (i.e., skip something like 'IMPRECISE' which has no =)
        # Refs:
        # http://stackoverflow.com/questions/4627981/creating-a-dictionary-from-a-string
        # http://stackoverflow.com/questions/17321138/python-one-line-list-comprehension-if-else-variants
        kv = [item.split('=') for item in info.split(';')]
        infodict = dict([v for v in kv if len(v)==2])

        chr2=infodict['CHR2']
        pos2=infodict['END']

        # We are only interested in translocations
        # if alt != "<TRA>": continue

        # Write the alphabetically smaller chr first
        if not options.both:
            if chr1 < chr2:
                o.write('\t'.join((chr1, pos1, chr2, pos2))+"\n")
            elif chr2 > chr1:
                o.write('\t'.join((chr2, pos2, chr1, pos1))+"\n")
            # if chr1 == chr2, write the smaller pos first
            elif pos2 > pos1:
                o.write('\t'.join((chr1, pos1, chr2, pos2))+"\n")
            else:
                o.write('\t'.join((chr2, pos2, chr1, pos1))+"\n")
        else:
            o.write('\t'.join((chr1, pos1, chr2, pos2))+"\n")
            o.write('\t'.join((chr2, pos2, chr1, pos1))+"\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] file.vaf ...
        """

    parser = OptionParser(usage_text, version="$Revision: 1.0 $")
    parser.add_option("-v", dest="verbose", action="store_true", help="verbose output")
    parser.add_option("-H", dest="header", action="store_true", help="Write out header")
    parser.add_option("-b", dest="both", action="store_true", help="Output pos1 and pos2")
    parser.add_option("-i", dest="infn", default="stdin", help="Input filename")
    parser.add_option("-o", dest="outfn", default="stdout", help="Output filename")

    (options, params) = parser.parse_args()
    if options.verbose:
        sys.stderr.write(str(options)+'\n')
        sys.stderr.write(str(params)+'\n')

    if options.infn == "stdin":
        f = sys.stdin 
    else:
        f = open(options.infn, 'r')
    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    if options.header:
        o.write('\t'.join("chr1", "pos1", "chr2", "pos2"))

    filterVAF(f, o, options)

    f.close()
    o.close()

if __name__ == '__main__':
    main()
