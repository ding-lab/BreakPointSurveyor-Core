# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript BreakpointCruncher.R [-v] [-P] [-A range.A] [-B range.B] [-F] [-g fn.ggp] [-p plot.type] 
#                [-a alpha][-c color][-f fill][-s shape][-z size] BP.fn breakpoint.ggp
#
#   Create or append various plots to breakpoint coordinate GGP file.  Chrom A coordinates are plotted
#   on X axis, B on Y.
#
#   Three types of plots supported:
#   * "point" - points are drawn at positions given by BPC coordinates
#   * "region" - rectangular region drawn at positions given by BPR data
#   * "segment" - linear segments drawn at positions given by BPR data
#
# Input arguments:
#
# * BP.fn is either BPC or BPR file, depending on argument of -p.
# * breakpoint.ggp is output file in GGP format
#
# Optional arguments:
#
# -g fn.ggp: Append graphics to given ggp file, rather than creating new ggp file.
# -A, -B: Define A and B range, resp., specified as "C" or "C:M-N", where C is chromosome name
#         and M,N are genomic start, stop positions.  Range is used for two distinct purposes:
#           1) Filter data before plotting
#           2) Specify plot region
#         Plot region is specified only if -g is not specified (i.e., creating new GGP)
#         Filtering data may be prevented with -F.
# -F: Do not filter data.  See above.
# -p: plot.type: one of "point", "region", "segment"; these require BPC, BPR, BPR data files, resp.  "point" is default.
# -P: Output as PDF file instead of GGP
#
# The following constant attributes can be defined:
#   -a alpha
#   -c color
#   -f fill
#   -s shape
#   -z size
# These correspond to attributes in ggplot geom call.  If -c or -f are not
# defined with a constant attribute, and an attribute column exists in the BPC
# or BPR file (see below), the value of color and fill (if appropriate) is
# passed as an aesthetic by the value of the attribute column.  (Note that need
# to define color_scale for specific colors, not implemented)
# see here for more details: http://docs.ggplot2.org/current/vignettes/ggplot2-specs.html
#
# Input data in TSV format, with the following columns: 
# * BPC: chromA, posA, chromB, posB, [attribute] 
# * BPR: chromA, posA.start, posA.end, chromB, posB.start, posB.end, [attribute] 
# The attribute column is optional.  See README.file_formats for more details

# Version 2.0: Abstracted out BPC and BPR plotting, appending to GGP

options("width"=180) # useful for debugging
library("bitops")
suppressPackageStartupMessages(library("ggplot2"))

# these control printing of genomic position labels
options(scipen=3)  # no scientific notation


# Return the command line argument associated with a given flag (i.e., -o foo),
# or the default value if argument not specified.
# Note that this will break if an argument is not supplied after the flag.
get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

# Return boolean specifying whether given flag appears in command line (i.e., -o),
get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

parse_chr = function(chrarg) {
    # Parse the chromosome limit string.  Accepted formats:
    # -c 14
    # -c 14:12345-12456
    range.pos <- NULL
    range.chr <- NULL
    if (chrarg != 'all') {
        chrlist<-strsplit(chrarg, ":")[[1]]
        range.chr <- strsplit(chrlist[1], ",")[[1]][1]  # chromosome name -- only one allowed
        if (length(chrlist) == 2) {
            range.pos<-as.numeric(strsplit(chrlist[2], "-")[[1]])
        }
    } else {
        print("Error: Must specify chromosome (-c).  Quitting")
        q()
    }

    return( list('range.pos'=range.pos, 'range.chr'=range.chr) )
}

# Usage: 
#   args = parse_args()
#   print(args$disease.filter)
parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    range.A = parse_chr(get_val_arg(args, "-A", "all"))  # accessible as range.pos.A, range.chr.A
    range.B = parse_chr(get_val_arg(args, "-B", "all"))  # accessible as range.pos.B, range.chr.B

    skip.data.filter = get_bool_arg(args, "-F")
    in.ggp = get_val_arg(args, "-g", NULL)
    plot.type = get_val_arg(args, "-p", "point")
    pdf.out = get_bool_arg(args, "-P")

    alpha = get_val_arg(args, "-a", NULL)
    color = get_val_arg(args, "-c", NULL)
    fill = get_val_arg(args, "-f", NULL)
    shape = get_val_arg(args, "-s", NULL)
    size = get_val_arg(args, "-z", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)];             args = args[-length(args)]
    data.fn = args[length(args)];      args = args[-length(args)]

    val = list('verbose'=verbose, 
            'range.pos.A'=range.A$range.pos, 'range.chr.A'=range.A$range.chr,  # note we are expanding out range.A
            'range.pos.B'=range.B$range.pos, 'range.chr.B'=range.B$range.chr, 
            'skip.data.filter'=skip.data.filter, 'in.ggp'=in.ggp, 'plot.type'=plot.type,
            'alpha'= alpha, 'color'= color, 'fill'= fill, 'shape'= shape, 'size'= size,
            'out.ggp'=out.ggp, 'data.fn' = data.fn, 'pdf.out'=pdf.out
            )
    if (val$verbose) { print(val) }

    return (val)
}

# reads in a BPC file which specifies breakpoint coordinates with optional attributes column.
read.BPC = function(BPC.fn) {
    data = read.table(BPC.fn, sep='\t')
    if (ncol(data) == 4) {
        names(data) = c("chrom.A", "pos.A", "chrom.B", "pos.B")
        data$attribute = "default"
    } else if (ncol(data) == 5) {
        names(data) = c("chrom.A", "pos.A", "chrom.B", "pos.B", "attribute")
    } else {
        stop("Unexpected number of rows in ", BPC.fn)
    }
    data$attribute = factor(data$attribute)
    return(data)
}

# * BPR: chrom.A, pos.A.start, pos.A.end, chrom.B, pos.B.start, pos.B.end, [attribute] 
read.BPR = function(BPR.fn) {
    data = read.table(BPR.fn, sep='\t')
    if (ncol(data) == 7) {
        names(data) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end")
        data$attribute = "default"
    } else if (ncol(data) == 8) {
        names(data) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "attribute")
    } else {
        stop("Unexpected number of rows in ", BPR.fn)
    }
    data$attribute = factor(data$attribute)
    return(data)
}

filter.BPC = function(data, range.chr.A, range.pos.A, range.chr.B, range.pos.B) {
    if (! is.null(range.chr.A)) data = data[data$chrom.A %in% range.chr.A,]
    if (! is.null(range.chr.B)) data = data[data$chrom.B %in% range.chr.B,]

    if (!is.null(range.pos.A)) {
        data = data[data$pos.A >= range.pos.A[1] & data$pos.A <= range.pos.A[2],]
    } 
    if (!is.null(range.pos.B)) {
        data = data[data$pos.B >= range.pos.B[1] & data$pos.B <= range.pos.B[2],]
    } 
    # Get rid of rows with NA in them.  It would be useful to refactor attributes as well.
    data = data[complete.cases(data),]
    return(data)
}

# evaluate whether there is any overlap between two segments a and b
# returns array of indices (which()).
is.overlapping = function(a.start, a.end, b.start, b.end) {
    # test whether a fully before b or b fully before a.  If neither, there is overlap.
    is.overlapping = which( ! (a.end < b.start | b.end < a.start))
    return(is.overlapping)
}

# We retain entire BPR region when any part of it falls within range of interest.  Region is
# not cropped.
# names(data) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end")
filter.BPR = function(data, range.chr.A, range.pos.A, range.chr.B, range.pos.B) {
    if (! is.null(range.chr.A)) data = data[data$chrom.A %in% range.chr.A,]
    if (! is.null(range.chr.B)) data = data[data$chrom.B %in% range.chr.B,]

    if (!is.null(range.pos.A)) {
        in.range.A = is.overlapping(data$pos.A.start, data$pos.A.end, range.pos.A[1], range.pos.A[2])
        data = data[in.range.A,]
    } 
    if (!is.null(range.pos.B)) {
        in.range.B = is.overlapping(data$pos.B.start, data$pos.B.end, range.pos.B[1], range.pos.B[2])
        data = data[in.range.B,]
    } 
    # Get rid of rows with NA in them.  It would be useful to refactor attributes as well.
    data = data[complete.cases(data),]
    return(data)
}

make.GGP = function(range.pos.A = NULL, range.pos.B = NULL) {
    # define basic properties of Breakpoint Coordinates GGP object.
    # TODO: see if this works.  It may be necessary to set some things at the end of plotting?
    ggp = ggplot() + xlab(NULL) + ylab(NULL) 
    
# going to have problems calling coord_cartesian multiple times here.  Better to build up argument
# as list and call coord_cartesian once.  There are several ways to go it...
# or do.call(runif,c(list(n=100),params)) 
# was:
# ggp = ggp + coord_cartesian(xlim=c(range.A.start, range.A.end), ylim=c(range.B.start, range.B.end))

    coord_cartesian.args = list()
    if (!is.null(range.pos.A)) coord_cartesian.args$xlim = c(range.pos.A[1], range.pos.A[2])
    if (!is.null(range.pos.B)) coord_cartesian.args$ylim = c(range.pos.B[1], range.pos.B[2])
    if (length(coord_cartesian.args) > 0) {
        ggp = ggp + do.call(coord_cartesian, coord_cartesian.args)
    }

    ggp = ggp + theme_bw()
    ggp = ggp + theme(axis.text=element_text(size=6), axis.text.y=element_text(angle=-90, hjust=0.5))  
#    ggp = ggp + theme(legend.position="none")  # is this what we want?

    # necessary?
    # ggp = ggp + scale_shape(guide=FALSE)  # get rid of shape legend but keep others
    return(ggp)
}


# default CTX values: geom_point.  color="#377EB8", alpha = 0.5
# rSBP: geom_point.  color="#4DAF4A", alpha=0.75, shape=3, size=4, show_guide=FALSE
# pindel region: geom_rect. fill="gray50", color=NA, alpha=0.2
#       when strand defined, aes(fill=strand)
# qSBP: geom_segment.  aes(color=contig.id), alpha=0.5

#   -a alpha
#   -c color
#   -f fill
#   -s shape
#   -z size
render.point = function(ggp, BPC, alpha=NULL, color=NULL, shape=NULL, size=NULL) {
# if alpha, shape, size, color are defined, their values are passed to corresponding arguments.  
# if color is not defined, then color is passed as an aesthetic with data$attribute as value.
# BPC: "chrom.A", "pos.A", "chrom.B", "pos.B", "attribute"


    # call to geom_point consists of three argument types:
    #   data = BPC
    #   aes arguments: mapping = aes( ... )
    #   static arguments, e.g., color = ..., 
    args = list()   # first collect all static arguments, then add aes and data as arguments
    aes.args = list(x="pos.A", y="pos.B")
    if (!is.null(alpha)) args$alpha = alpha
    if (!is.null(shape)) args$shape = shape
    if (!is.null(size)) args$size = size
    if (!is.null(color)) args$color = color
    else { # if not specified, color is an aes with value given by attribute column
        aes.args$color = "attribute"
    }
    args$data=BPC

    # construct aes argument from strings
    # http://tolstoy.newcastle.edu.au/R/e3/help/07/12/6372.html
    args$mapping = do.call(aes_string, aes.args)  

    # now call geom_point with all the arguments
    ggp = ggp + do.call(geom_point, args)  
    return(ggp)
}

# Render region based on BPR data
# names(BPR) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "attribute")
render.region = function(ggp, BPR, alpha=NA, color=NA, fill=NA) {
    if (! is.null(BPR) && nrow(BPR)>0) {
        # if fill == NA, aes(fill=attribute)
        # if color == NA, aes(color=attribute)
        # if both true, aes(fill=attribute, color=attribute)
        # define argument to geom_rect fill = fill if fill != NA
        # define argument to geom_rect color = color if color != NA
        ggp = ggp + geom_rect(data=BPR, mapping=aes(xmin=pos.A.start, xmax=pos.A.end, ymin=pos.B.start, ymax=pos.B.end), fill=fill, color=color, alpha=alpha)
    }

    return(ggp)
}

render.segment = function(ggp, BPR, color=NA, alpha=NA, size=NA, linetype=NA) {  
    if (! is.null(BPR) && nrow(BPR)>0) {
        qSBP.data$contig.id = factor(qSBP.data$contig.id)
        p = p + geom_segment(data = qSBP.data, aes(x=h.bp.rpos.bpA, y=v.bp.rpos.bpA, xend=h.bp.rpos.bpB, yend=v.bp.rpos.bpB, color=contig.id ), alpha=0.5)
    }
}

args = parse_args()

if (is.null(args$in.ggp)) {
    ggp = make.GGP(args$range.pos.A, args$range.pos.B)
} else {
    ggp = readRDS(args$ggp.in)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
}

if (args$plot.type == "point") {
    breakpoints.bpc = read.BPC(args$data.fn)

    if (!args$skip.data.filter) {
        breakpoints.bpc = filter.BPC(breakpoints.bpc, args$range.chr.A, args$range.pos.A, args$range.chr.B, args$range.pos.B)
    }
    ggp = render.point(ggp, breakpoints.bpc, alpha=args$alpha, color=args$color, shape=args$shape, size=args$size)

} else if (args$plot.type == "region") {
#    filter.BPR(data, range.chr.A, range.pos.A, range.chr.B, range.pos.B) 
    print("Unimplemented")

} else if (args$plot.type == "segment") {
    print("Unimplemented")

} else {
    stop("Unknown plot type", args$plot.type)
}

# TODO:
# test and develop:
#   * adding layers to a GGP object
#   * Pindel regions
#   * segments

# now save breakpoints.ggp to file
# http://stat.ethz.ch/R-manual/R-devel/library/base/html/save.html
if (!args$pdf.out) {
    cat(paste("Saving to GPP file", args$out.ggp, "\n"))
    saveRDS(ggp, file=args$out.ggp)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
} else {
    cat(paste("Saving to PDF file", args$out.ggp, "\n"))
    ggsave(plot=ggp, filename=args$out.ggp, useDingbats=FALSE)
}
