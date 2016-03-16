#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu

# Cluster Breakpoints into Breakpoint regions.

# BPC file is read with breakpoint coordinates.
# For a given chromosome pair, all breakpoints within a given distance of each other
# are collected into individual clusters.  Such clusters then define the regions of
# interest as written to a BPR file.
# Implementation details illustrated in BreakpointSurveyor/doc/BPS_clustering.pdf

# Details on data formats:
#     In general, we represent breakpoints as a coordinate given by a pair of chrom/pos (Breakpoint Coordinates, or BPC).
#     Alternatively, we may consider a breakpoint region (BPR), which has a pair of chrom/pos.start/pos.end values.
#     Breakpoints with precise positions (e.g., discordant pair positions) will be represented by the former,
#     while regions such as SV Events will be represented by the latter.
# 
#     Each breakpoint coordinate or region is represented just once, with chromA < chromB, or posA < posB if chromA==chromB
# 
#     Both are represented in a TSV file,
#       BPC: chromA, posA, chromB, posB, [attribute]
#       BPR: chromA, posA.start, posA.end, chromB, posB.start, posB.end, [attribute]
#     Both datatypes have optional attribute column

import sys

class Breakpoints:
    """Read, filter, and store breakpoint data from BPC file."""
    # internally, breakpoints are represented as tuples, (chromA, posA, chromB, posB)
    # where chromA, chromB are the chromosomes of interest as defined by user
    # the index into list of breakpoints is the bp_id and uniquely identifies a breakpoint

    def __init__(self, f, chromA, chromB, noSwap):
        self.breakpoints = []
        for line in f:
            if line.startswith("#"): continue
            chrom_1, pos_1, chrom_2, pos_2 = line.rstrip().split("\t")
            if chrom_1 == chromA and chrom_2 == chromB:
                self.breakpoints.append( (chrom_1, int(pos_1), chrom_2, int(pos_2)) )
            elif chrom_2 == chromA and chrom_1 == chromB and not noSwap:
                self.breakpoints.append( (chrom_2, int(pos_2), chrom_1, int(pos_1)) )

    def __repr__(self):
        return str(self.breakpoints)

    def overlapping(self, bp1_id, bp2_id, R):
        """Test whether breakpoint 1 is within distance radius of breakpoint 2 in any direction"""
        # algorithm based partly on that here: http://stackoverflow.com/questions/13390333/two-rectangles-intersection
        # bp1 has coordinates (A1, B1), bp2 has (A2, B2)
        bp1 = self.breakpoints[bp1_id]
        bp2 = self.breakpoints[bp2_id]
        A1, B1 = bp1[1], bp1[3]
        A2, B2 = bp2[1], bp2[3]
        overlap = True
        if A1 + R < A2 - R: overlap = False
        if A1 - R > A2 + R: overlap = False
        if B1 + R < B2 - R: overlap = False
        if B1 - R > B2 + R: overlap = False
        return overlap

    def getBPList(self): 
        return self.breakpoints

    def BP(self, i): 
        return self.breakpoints[i]

    def getBPCount(self): 
        return len(self.breakpoints)

    def getBPRange(self, c):
        "Given list of breakpoints ids c, return range start and stop for chrom A and B"
        startA, startB, endA, endB, chromA, chromB = float("inf"), float("inf"), float("-inf"), float("-inf"), None, None
        for bp in c:
            (chromA, posA, chromB, posB) = self.breakpoints[bp]
            startA, startB = min(startA, posA), min(startB, posB)
            endA, endB = max(endA, posA), max(endB, posB)
        return chromA, startA, endA, chromB, startB, endB


class Clusters:
    """Maintains, examines, and manipulates a list of clusters, each of which holds overlapping breakpoints."""
    def __init__(self):
        # the principal data structure is a dictionary of lists which will hold breakpoint IDs.
        self.clusters = {}
        # cluster IDs are incrementing integers
        self.next_cid = 0

    def __repr__(self):
        return str(self.clusters)

    def newCluster(self, bpid=None):
        "Create new cluster, optionally add bpid, and return its cluster id"
        new_cid = self.next_cid
        self.next_cid += 1
        self.clusters[new_cid] = []
        if bpid is not None: self.clusters[new_cid].append(bpid)
        return new_cid

    def add(self, cid, bpid):
        "Add breakpoint id to given cluster id"
        self.clusters[cid].append(bpid)

    def getCid(self, bpid):
        "Return cluster ID containing given bp id, or None if bp is not in any cluster"
        for cid in self.clusters.keys():
            if bpid in self.clusters[cid]: return cid
        return None

    def merge(self, cidA, cidB):
        "Merge clusters A and B.  All bp in cluster B added to cluster A, and cluster B is removed"
        self.clusters[cidA] += self.clusters[cidB]
        del self.clusters[cidB]


    def writeBPR(self, o, bp, header=False, countBP=False):
        "Write contents of this object in BPR format"
        # BPR object has columns, chromA, posA.start, posA.end, chromB, posB.start, posB.end
        if header: 
            if countBP:
                o.write('\t'.join( ("chromA", "startA", "endA", "chromB", "startB", "endB", "breakpointCount") )+"\n")
            else:
                o.write('\t'.join( ("chromA", "startA", "endA", "chromB", "startB", "endB") )+"\n")
        for k,c in self.clusters.iteritems():  # c is a list of breakpoint ids
            if countBP:
                o.write('\t'.join( map(str, bp.getBPRange(c) + (len(c),) ) )+"\n")
            else:
                o.write('\t'.join( map(str, bp.getBPRange(c)) )+"\n")

def makeBreakpointClusters(breakpoints, radius):
    """Merge into clusters those breakpoints which are within radius distance of each other.  Return clusters object."""
    # see notes BreakpointSurveyor/doc/BPS_clustering.pdf for algorithm details
    clusters = Clusters()
    for bp_i in range(0, breakpoints.getBPCount()):  # bp_i loops over all breakpoints
        currentCluster = clusters.getCid(bp_i)
        if currentCluster is None:
            currentCluster = clusters.newCluster(bp_i)
        for bp_j in range(bp_i+1, breakpoints.getBPCount()):  # bp_j loops over i..N breakpoints
            if breakpoints.overlapping(bp_i, bp_j, radius):
                cluster_j = clusters.getCid(bp_j)
                if cluster_j is not None and cluster_j is not currentCluster:
                    clusters.merge(currentCluster, cluster_j)
                else:
                    if cluster_j is not currentCluster:
                        clusters.add(currentCluster, bp_j)
    return clusters
        


def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] BPC.dat BPR.dat ...
Cluster Breakpoints into Breakpoint regions.

BPC file is read with breakpoint coordinates.  For a given chromosome pair, all
breakpoints within a given distance of each other along both chromosomes are
collected into collective clusters.  Such clusters then define the regions of
interest as written to a BPR file."""

    parser = OptionParser(usage_text, version="$Revision: 1.2 $")
    parser.add_option("-R", dest="radius", default="10000000", help="Cluster all breakpoints within this distance of each other.  Default 10M.")
    parser.add_option("-A", dest="chromA", default=None, help="First chromosome of interest")
    parser.add_option("-B", dest="chromB", default=None, help="Second chromosome of interest")
    parser.add_option("-n", dest="noSwap", action="store_true", help="Only evaluate (chromA,chromB) and ignore (chromB,chromA)")
    parser.add_option("-H", dest="header", action="store_true", help="Print BPR header")
    parser.add_option("-d", dest="debug", action="store_true", help="Print cluster details as comments in output file")
    parser.add_option("-c", dest="countBP", action="store_true", help="Print count of breakpoints per cluster as final BPR column")

    (options, params) = parser.parse_args()

    if (len(params) != 2):
        parser.error("Pass two arguments.")
    infn = params[0]
    outfn = params[1]

    if infn == "stdin":
        f = sys.stdin
    else:
        f = open(infn, 'r')
    if outfn == "stdout":
        o = sys.stdout
    else:
        o = open(outfn, "w")

    breakpoints = Breakpoints(f, options.chromA, options.chromB, options.noSwap) 
    clusters = makeBreakpointClusters(breakpoints, int(options.radius))

    if (options.debug):
        o.write( "# Number breakpoints: " + str(breakpoints.getBPCount()) + "\n")
        o.write( "# Cluster radius: " + options.radius + "\n" )
        o.write( "# Cluster dump: "+ str(clusters) + "\n" )

    clusters.writeBPR(o, breakpoints, options.header, options.countBP)
    f.close()
    o.close()


if __name__ == '__main__':
    main()
