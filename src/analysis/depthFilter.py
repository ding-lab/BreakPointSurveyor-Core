#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
#
# depthFilter.py: Evaluate read depth in subsection of chrom, subsampled to give reasonable output size, 
#                 optimized for performance
#
# V3.2.  Discover chrom length automatically
# V3.1.  Option to get depth by region or basepair implemented
# V3.0.  Rewritten to use pysam 
# V2.1.  Simple steps to optimize this.  More tips here: https://www.python.org/doc/essays/list2str/
# V2.0.  Dropped list support (see Q_ReadDepth/bin/depthFilter.py).  Support for multiple segments
#        per BED file.  Stride automatically determined for each.
# V1.3.  Bug fix (stride badly implemented)
# V1.2.  Support for filtering by specific positions
# V1.1.  Added support for automatic chr length detection.

# from http://pysam.readthedocs.org/en/latest/api.html
# Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.
# we retain 0-based coordinates in the position output.  Note that pysam count_coverage differs on occasion from 
# samtools depth by one or two reads at a given position.  They don't seem to be precisely the same, but this
# difference does not seem significant and will not be pursued.


# Performance and pysam
# The typical approach to getting read depth is to evaluate it for an entire region and then subsample
# to obtain the "npts" number of points.
# Preliminary testing indicated that this by-region approach for 1000 points takes about 10 sec/1Mbp region.
# however, large regions (tens of Mbp) cause out of memory errors.
# an alternative approach for large regions is to evaluate depth at a single position (per-basepair approach).
# this reduces the memory required at the expense of more calls to read depth (increasing I/O).  To strike
# a balance between these two approaches, we use per-region analysis to get read depth for small regions, and 
# a per-basepair approach to get read depth for large regions.  
# We evaluate stride as region_size / npts, where region_size = pos_end - pos_start, and npts is the number
# of samples requested.  We then use a stride_threshold to determine which approach to take:
# for stride > stride_threshold, use per-region depth analysis; otherwise, use per_basepair analysis.

import sys, os
import pysam
import numpy
import timeit  # http://stackoverflow.com/questions/15707056/get-time-of-execution-of-a-block-of-code-in-python-2-7


def openAlignment(fn):
    """
    Opens SAM, BAM, or CRAM file and returns pysam.AlignmentFile object.
    Filetype is detected based on extension.
    """
    # This can be useful for understanding how things work:
    #   https://github.com/pysam-developers/pysam/blob/master/pysam/calignmentfile.pyx#L450
    filename, file_extension = os.path.splitext(fn)
    if file_extension == ".sam":
        alignmentFile = pysam.AlignmentFile( fn, "r" )  
    elif file_extension == ".bam":
        alignmentFile = pysam.AlignmentFile( fn, "rb" )
    elif file_extension == ".cram":
        alignmentFile = pysam.AlignmentFile( fn, "rc" )
    else:
        raise Exception("Unknown input filename extension "+file_extension)
    return alignmentFile

def getStride(start, end, npts):
    if npts == 0:  # npts == 0 is a special case, so that every line is printed
        s = 1
    else:
        # obtain stride so about npts lines are printed per segment
        # since stride needs to be an integer, estimate closest stride and the number of points it gives.
        s = int(round(float(end - start)/float(npts)))
        if s == 0: s = 1
    return s

def getDepthByRegion(alignmentFile, chrom, start, end, stride):
    """
    Evaluate the read depth in alignmentFile between start and end positions, subsampled every stride postions.
    Return list depth
    """
    depthACGT = alignmentFile.count_coverage(reference=chrom, start=start, end=end)
    depth = sum(map(numpy.array, depthACGT))[::stride]
    return depth

def getDepthByBasepair(alignmentFile, chrom, pos):
    """
    Evaluate read depth at each position as given by list pos and chrom.
    return depth at those positions.
    """
    depth = numpy.array(pos)
    for i, p in enumerate(pos):
        depthACGT = alignmentFile.count_coverage(reference=chrom, start=p, end=p+1)
        depth[i] = sum(map(numpy.array, depthACGT))
    return depth

def getDepth(alignmentFile, chrom, start, end, npts, stride_threshold):
    """
    Evaluate the read depth in alignmentFile between start and end positions, subsampled every stride postions.
    If byRegion, evaluate depth in entire genomic region first and subsample second; otherwise
    evaluate depth iteratively at each position of interest.

    Return tuple (pos, depth), with pos and depth arrays of the same length indicating chrom position and depth.
    """
    stride = getStride(start, end, npts)
    pos = numpy.arange(start, end, stride)
    print "Stride:", stride
    if stride < stride_threshold:
        print "By Region"
        start_time = timeit.default_timer()
        depth = getDepthByRegion(alignmentFile, chrom, start, end, stride)
        print timeit.default_timer() - start_time, "sec"
    else:
        print "By Basepair"
        start_time = timeit.default_timer()
        depth = getDepthByBasepair(alignmentFile, chrom, pos)
        print timeit.default_timer() - start_time, "sec"
    return (pos, depth)


def writeDepth(chrom, pos, depth, o):
    # sample 'samtools depth' line:
    #   14  68649199    5
    for p, d in zip(pos, depth):
        o.write( "\t".join( (chrom, str(p), str(d)) )+"\n")

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] chrom start end fn ...
        Parse alignment file and output read depth in format, "chrom pos depth"
        chrom: name of reference sequence.  
        start, end: start and end position.  end="END" obtains end position from reference
        fn: alignment filename. Accepted fn file formats (sam, bam, cram) determined by filename extension.
        """
    parser = OptionParser(usage_text, version="$Revision: 3.0 $")
    parser.add_option("-v", dest="verbose", action="store_true", help="verbose output")
    parser.add_option("-N", dest="npts", default="0", help="Subsample depth to output depth at approximately npts equally spaced postions.")
    parser.add_option("-o", dest="outfn", default="stdout", help="output filename")
    parser.add_option("-t", dest="timing", action="store_true", help="Output timing info to stderr")
    parser.add_option("-S", dest="stride_threshold", default=500, help="Evaluate depth by region, by position if stride less, greater than this value, resp.")

    (options, params) = parser.parse_args()
    if options.verbose:
        sys.stderr.write(str(options)+'\n')
        sys.stderr.write(str(params)+'\n')
    if (len(params) != 4):
        parser.error("Please pass BED and output filenames.")

    chrom, start, end, fn = params[0], int(params[1]), params[2], params[3]

    if options.outfn == "stdout":
        o = sys.stdout
    else:
        o = open(options.outfn, "w")

    print "opening", fn
    f = openAlignment(fn)

    if chrom not in f.references:
        parser.error("Unknown chrom name.  List of known chrom:\n" + str(f.references))

    # discover chrom length if necessary
    if end == "END":
        end = f.lengths[f.references.index(chrom)]
    else:
        end = int(end)

    start_time = timeit.default_timer()
    (pos, depth) = getDepth(f, chrom, start, end, int(options.npts), int(options.stride_threshold))
    elapsed = timeit.default_timer() - start_time

    writeDepth(chrom, pos, depth, o)
    if options.timing:
        sys.stderr.write("Elapsed time %0.3f sec\n" % elapsed)

if __name__ == '__main__':
    main()
