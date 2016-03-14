# DNAcopy necessary for CBS
library(DNAcopy)  # DNAcopy is installed with source("http://bioconductor.org/biocLite.R"); biocLite(); biocLite("DNAcopy")

normalizeCopyNumber = function(depth.list, num.reads, read.length) {
    # Based on CombinedPlot/CombinedPlot.old/orig2013/src/depth_plot.R
    # given read depth data, normalize by average read depth to obtain copy number
    # num_reads is the number of reads in BAM (either total or mapped, your call), typically obtained from .flagstat file
    # read_length is the average read length, also obtained from BAM
    # we normalize the read depth by the average read depth across the genome, defined as 2 * num_reads * read_length / genome_length
    #  - factor 2 is so that copy number comes out to be 2 on average
    num.bp.genome = 3137144693.
    depth.factor = num.reads * read.length / num.bp.genome
#    print(sprintf("Average BAM read depth = %f", depth.factor))
    norm.depth = 2. * depth.list / depth.factor
    return(norm.depth)
}

# return normalized read length based on values of doLog, num.reads, and read.length
# - if doLog=TRUE, norm.depth = log2(depth/<depth>), where depth is read depth and <depth> is the median depth of the filtered dataset
# - if doLog=FALSE:
#   * if num.reads and read.length are specified, norm.depth is read depth divided by average read depth across the genome (as obtained from number of reads and read length)
#   * if num.reads or read.length not specified, norm.depth is depth unchanged.
# return list ['norm.depth', 'method.name']
normalize.depth = function(depth.list, doLog, num.reads=NA, read.length=NA) {
    if (doLog) {
        # evaluate log.depth as log2(depth/<depth>)
        norm.depth = log2(depth.list/median(depth.list,na.rm=T))
        method="Log Depth Ratio"
    } else if (!is.na(num.reads) & (!is.na(read.length))) {
        norm.depth = normalizeCopyNumber(depth.list, num.reads, read.length)
        method="Copy Number"
    } else {
        norm.depth=depth.list
        method="Read Depth" 
    }
    return(list('norm.depth'=norm.depth, 'method'=method))
}

read.depth = function(depth.fn) {
    depth = read.table(depth.fn, header=FALSE, sep="\t", colClasses=c("character","numeric","numeric"), col.names=c("chrom", "pos", "depth"))
    return(depth)
}

filter.depth = function(depth, range.chr, range.pos) {
    if (! is.null(range.chr)) depth = depth[depth$chrom %in% range.chr,]

    if (!is.null(range.pos)) {
        depth = depth[depth$pos >= range.pos[1] & depth$pos <= range.pos[2],]
    } 
    # Get rid of rows with NA in them.  It would be useful to refactor attributes as well.
    depth = depth[complete.cases(depth),]
    return(depth)
}

calculateCBS = function(depth, do.log) {
    #CNA.object = CNA( log.depth = log2(depth$depth/median(depth$depth,na.rm=T)), chrom = depth$chrom, maploc = depth$pos, data.type = 'logratio')

    # Documentation is not clear whether non-logratio data supported.
    #   http://bioconductor.org/packages/release/bioc/html/DNAcopy.html
    # If we are confident it is, then should process norm.depth so that copy number scaling also supported. 
    if (do.log) {
        CNA.object = CNA( genomdat = log2(depth$depth/median(depth$depth,na.rm=T)), chrom = depth$chrom, maploc = depth$pos)
    } else {
        CNA.object = CNA( genomdat = depth$norm.depth, chrom = depth$chrom, maploc = depth$pos)
    }
    CNA.object = smooth.CNA(CNA.object)  
    d = segment(CNA.object, verbose=0, min.width=5, undo.splits="sdundo", undo.SD=3, alpha=1e-5)  # Smaller alpha makes fewer segments
    segs = d$output
    #write.table(segs[,2:6], file=outfn, row.names=F, col.names=F, quote=F, sep="\t")
    return (segs[,2:6])
}

# Construct CBS data frame which segments read depth into discrete levels
# Underlying code seems to require log2 scaling, although there is experimental support for segmentation of raw read counts.
getCBS = function(depth, do.log) {
    # Construct CBS data frame which segments read depth into discrete levels
    CBS = calculateCBS(depth, do.log)
    names(CBS)[names(CBS)=="loc.start"] = 'start' 
    names(CBS)[names(CBS)=="loc.end"] = 'end'
#      chrom loc.start  loc.end num.mark seg.mean g data.type
#      1    14  68649196 68654404      373  -0.5861       1       CBS     # outdated

    CBS$num.mark = NULL
    CBS$start = as.numeric(CBS$start)
    CBS$end = as.numeric(CBS$end)
    CBS$pos = round(CBS$start + (CBS$end-CBS$start)/2)
    CBS$g = rownames(CBS)  # This becomes the group used for geom_segment, so they are disconnected.
    return(CBS)
}
