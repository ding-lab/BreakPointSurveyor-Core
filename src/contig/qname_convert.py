#!/usr/bin/python
# Matthew A. Wyczalkowski, m.wyczalkowski@wustl.edu

# Usage:
# cat long.sam | python qname_convert.py > md5.sam 

import hashlib, sys

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] ...
        Convert QNAME in SAM file to MD5 hash.  This may be necessary in cases where the name of a contig is too long for SAMTOOLS to parse.
        """

    parser = OptionParser(usage_text, version="$Revision: 1.0 $")

    (options, params) = parser.parse_args()

    f = sys.stdin

    for line in f:
        if line[0] == "@": 
            print line.rstrip()
            continue

        t = line.rstrip().split("\t")
        t[0] = hashlib.md5(t[0]).hexdigest()
        print '\t'.join(t)


if __name__ == '__main__':
    main()
