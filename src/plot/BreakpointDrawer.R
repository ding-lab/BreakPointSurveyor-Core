# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript BreakpointDrawer.R [-v] [-P] [-l] [-A range.A] [-B range.B] [-F] [-G fn.ggp] [-p plot.type] [-M]
#                [-a alpha][-c color][-f fill][-s shape][-t linetype][-z size] BP.fn breakpoint.ggp
#
#   Create or append various plots to breakpoint coordinate GGP file.  Chrom A coordinates are plotted
#   on X axis, B on Y.
#
#   Three types of plots supported:
#   * "point" - points are drawn at positions given by BPC coordinates
#   * "region" - rectangular region drawn at positions given by BPR data
#   * "segment" - linear segments drawn at positions given by BPR data.  Unimplemented
#
# Input arguments:
#
# * BP.fn is either BPC or BPR file, depending on argument of -p.
# * breakpoint.ggp is output file in GGP format. If -P defined, .pdf appended to filename if necessary
#
# Optional arguments:
#
# -G fn.ggp: Append graphics to given ggp file, rather than creating new ggp file.
# -A, -B: Define A and B range, resp., specified as "C" or "C:M-N", where C is chromosome name
#         and M,N are genomic start, stop positions.  Range is used for two distinct purposes:
#           1) Filter data before plotting
#           2) Specify plot region
#         Plot region is specified only if -G is not specified (i.e., creating new GGP)
#         Filtering data may be prevented with -F.
# -F: Do not filter data.  See above.
# -l: Swap A, B columns when reading in BPC or BPR, so chrom A > chrom B by string comparison
# -p: plot.type: one of "point", "region", "segment"; these require BPC, BPR, BPR data files, resp.  "point" is default.
# -P: Output as PDF file instead of GGP.  This is primarily for convenience and debugging.
# -M: format genomic coordinates without commas 
#
# The following constant attributes can be defined:
#   -a alpha
#   -c color
#   -f fill
#   -s shape
#   -t linetype
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
library(scales)    # necessary for commas
suppressPackageStartupMessages(library("ggplot2"))

# these control printing of genomic position labels
options(scipen=3)  # no scientific notation

source_relative = function(source.fn) {
    # http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
    initial.options = commandArgs(trailingOnly = FALSE)
    file.arg.name = "--file="
    script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename = dirname(script.name)
    other.name = paste(sep="/", script.basename, source.fn)
#    print(paste("Sourcing",other.name,"from",script.name))
    source(other.name)
}
source_relative("BPS_PlotUtil.R")
source_relative("../util/BPS_Util.R")

# Usage: 
#   args = parse_args()
#   print(args$disease.filter)
parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    range.A = parse.range.str(get_val_arg(args, "-A", NULL))  # accessible as range.pos.A, range.chr.A
    range.B = parse.range.str(get_val_arg(args, "-B", NULL))  # accessible as range.pos.B, range.chr.B

    skip.data.filter = get_bool_arg(args, "-F")
    in.ggp = get_val_arg(args, "-G", NULL)
    plot.type = get_val_arg(args, "-p", "point")
    pdf.out = get_bool_arg(args, "-P")
    no.commas = get_bool_arg(args, "-M")
    flip.ab = get_bool_arg(args, "-l")

    alpha = as.numeric(get_val_arg(args, "-a", NA))
    shape = as.numeric(get_val_arg(args, "-s", NA))
    size = as.numeric(get_val_arg(args, "-z", NA))
    color = get_val_arg(args, "-c", NULL)
    fill = get_val_arg(args, "-f", NULL)
    linetype = get_val_arg(args, "-t", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)];             args = args[-length(args)]
    data.fn = args[length(args)];      args = args[-length(args)]

    val = list('verbose'=verbose, 
            'range.pos.A'=range.A$range.pos, 'range.chr.A'=range.A$range.chr,  # note we are expanding out range.A
            'range.pos.B'=range.B$range.pos, 'range.chr.B'=range.B$range.chr, 
            'skip.data.filter'=skip.data.filter, 'in.ggp'=in.ggp, 'plot.type'=plot.type,
            'alpha'= alpha, 'color'= color, 'fill'= fill, 'shape'= shape, 'size'= size, 'linetype'=linetype,
            'out.ggp'=out.ggp, 'data.fn' = data.fn, 'pdf.out'=pdf.out, 'no.commas'=no.commas, 'flip.ab'=flip.ab
            )
    if (val$verbose) { print(val) }

    return (val)
}

make.breakpoint.GGP = function(range.pos.A = NULL, range.pos.B = NULL, no.commas=FALSE) {
    # define basic properties of Breakpoint Coordinates GGP object.
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

    labels = if (no.commas) waiver() else comma
    ggp = ggp + scale_y_continuous(labels=labels, trans=reverse_trans())  # make y position increase in the downward direction
    ggp = ggp + scale_x_continuous(labels=labels)  

    return(ggp)
}

# for defining colors... based on past results.
# default CTX values: geom_point.  color="#377EB8", alpha = 0.5
# rSBP: geom_point.  color="#4DAF4A", alpha=0.75, shape=3, size=4, show_guide=FALSE
# pindel region: geom_rect. fill="gray50", color=NA, alpha=0.2
#       when strand defined, aes(fill=strand)
# qSBP: geom_segment.  aes(color=contig.id), alpha=0.5

render.point = function(ggp, BPC, alpha=NA, color=NULL, fill=NULL, shape=NA, size=NA) {
# if alpha, shape, size, color are defined, their values are passed to corresponding arguments.  
# if color is not defined, then color is passed as an aesthetic with data$attribute as value.
# BPC: "chrom.A", "pos.A", "chrom.B", "pos.B", "attribute"

    # call to geom_point consists of three argument types:
    #   data = BPC
    #   aes arguments: mapping = aes( ... )
    #   static arguments, e.g., color = ..., 
    args = list()   # first collect all static arguments, then add aes and data as arguments
    aes.args = list(x="pos.A", y="pos.B")
    if (!is.na(alpha)) args$alpha = alpha
    if (!is.na(shape)) args$shape = shape
    if (!is.na(size) != 0) args$size = size
    if (!is.null(color)) args$color = color
    if (!is.null(fill)) args$fill = fill
    else { # if not specified, color is an aes with value given by attribute column
        aes.args$color = "attribute"
    }
    args$data=BPC

    # construct aes argument from strings
    # http://tolstoy.newcastle.edu.au/R/e3/help/07/12/6372.html
    args$mapping = do.call(aes_string, aes.args)  

    # now call geom_point with all the arguments
    ggp = ggp + do.call(geom_point, args)  

    # If there was no attribute column specified in BPC file (all attributes are "default") then
    # get rid of color legend.
    # TODO: Getting rid of legend if it is default data is good to do, but it will clobber any existing legends so
    # seems more appropriate to do later, or have that be defined by a flag.
#    attrib.uniq = unique(BPC$attribute)
#    if (length(attrib.uniq) == 1 & attrib.uniq[1] == "default") {
#        ggp = ggp + scale_color_discrete(guide=FALSE)  
#    }

    return(ggp)
}

# Render region based on BPR data
# names(BPR) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "attribute")
render.region = function(ggp, BPR, alpha=NA, color=NULL, fill=NULL, size=NA) {
    # if alpha, fill, or color are defined, their values are passed to corresponding arguments.  
    # if color is not defined, then color is passed as an aesthetic with data$attribute as value.
    if (is.null(BPR)) return (ggp)

    args = list()   # first collect all static arguments, then add aes and data as arguments
    aes.args = list(xmin="pos.A.start", xmax="pos.A.end", ymin="pos.B.start", ymax="pos.B.end")

    if (!is.na(alpha)) args$alpha = alpha
    if (!is.na(size)) args$size = size
    if (!is.null(fill)) args$fill = fill
    if (!is.null(color)) args$color = color
    else { # if not specified, color is an aes with value given by attribute column
        aes.args$color = "attribute"
    }
    args$data=BPR
    args$mapping = do.call(aes_string, aes.args)  
    ggp = ggp + do.call(geom_rect, args)  

    # If there was no attribute column specified in BPR file (all attributes are "default") then
    # get rid of color legend.
    attrib.uniq = unique(BPR$attribute)
    if (length(attrib.uniq) == 1 & attrib.uniq[1] == "default") {
        ggp = ggp + scale_color_discrete(guide=FALSE)  
    }

    # basic call: geom_rect(data=BPR, mapping=aes(xmin=pos.A.start, xmax=pos.A.end, ymin=pos.B.start, ymax=pos.B.end), fill=fill, color=color, alpha=alpha)
    return(ggp)
}

render.segment = function(ggp, BPR, color=NULL, alpha=NA, size=NA, linetype=NULL) {  
    # TODO: Implement color, alpha, size, linetype selection
    if (! is.null(BPR) && nrow(BPR)>0) {
        qSBP.data$contig.id = factor(qSBP.data$contig.id)
        p = p + geom_segment(data = qSBP.data, aes(x=A.bp.rpos.bpM, y=B.bp.rpos.bpM, xend=A.bp.rpos.bpN, yend=B.bp.rpos.bpN, color=contig.id ), alpha=0.5)
    }
}

args = parse_args()

if (is.null(args$in.ggp)) {
    ggp = make.breakpoint.GGP(args$range.pos.A, args$range.pos.B, args$no.commas)
} else {
    ggp = readRDS(args$in.ggp)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
}

if (args$plot.type == "point") {
    breakpoints.bpc = read.BPC(args$data.fn, args$flip.ab)

    if (!args$skip.data.filter) {
        breakpoints.bpc = filter.BPC(breakpoints.bpc, args$range.chr.A, args$range.pos.A, args$range.chr.B, args$range.pos.B)
    }
    cat(paste("Plotting ", nrow(breakpoints.bpc), " points\n"))
    ggp = render.point(ggp, breakpoints.bpc, alpha=args$alpha, color=args$color, fill=args$fill, shape=args$shape, size=args$size)

} else if (args$plot.type == "region") {
    breakpoints.bpr = read.BPR(args$data.fn, args$flip.ab)

    if (!args$skip.data.filter) {
        breakpoints.bpr = filter.BPR(breakpoints.bpr, args$range.chr.A, args$range.pos.A, args$range.chr.B, args$range.pos.B)
    }
    cat(paste("Plotting ", nrow(breakpoints.bpr), " regions\n"))
    ggp = render.region(ggp, breakpoints.bpr, alpha=args$alpha, color=args$color, fill=args$fill, size=args$size)

} else if (args$plot.type == "segment") {
    print("Unimplemented")

} else {
    stop("Unknown plot type", args$plot.type)
}

write.GGP(ggp, args$out.ggp, args$pdf.out)
