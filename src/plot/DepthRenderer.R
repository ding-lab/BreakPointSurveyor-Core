# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript DepthRenderer.R [-v] [-P] [-A range] [-F] [-G fn.ggp] [-p plot.type] 
#                [-u num.reads] [-l read.length] [-m chrom] [-L] [-B]
#                [-a alpha] [-c color] [-f fill] [-s shape] [-z size] [-t linetype] data.fn out.ggp
#
#   These values of plot.type are supported:
#   * "depth" - read depth is plotted as a sequence of points.  "depth" file format
#   * "CBS" - Plot segments of uniform depth (circular binary segmentation). "depth" file format
#   * "point" - points are drawn at pos given by posA. y-position is jittered, posB ignored. BPC input
#   * "vline" - vertical lines drawn at pos given by posA.  posB ignored.  BPC input
#   * "region" - region at posA.start, posA.end with indefinite Y extent.  BPR input
#   * "segment" - disconnected segments with X, Y endpoints given by A, B BPR regions, resp. (Used for precalculated CBS) 
#
# Input arguments:
#
# * data.fn is depth, BPC, or BPR file, depending on argument of -p.
# * out.ggp is filename of GGP output  If -P defined, .pdf appended to filename if necessary
#
# Optional arguments:
#
# -L: Plot read depth as log2(depth / <depth>), where <depth> is the median read depth
# -u: number of reads (typically mapped) in BAM file.  Necessary for normalizing to copy number
# -l: read length.  Necessary for normalizing to copy number
#
# -m chrom: if plotting point or region from BPC or BPR data, resp., define whether using chrom A or B.  "A" is default.
# -A: Define range, specified as "C" or "C:M-N", where C is chromosome name
#         and M,N are genomic start, stop positions.  Range is used for two distinct purposes:
#           1) Filter data before plotting
#           2) Specify plot region
#         Plot region is specified only if -G is not specified (i.e., creating new GGP)
#         Filtering data may be prevented with -F.  Applies to chrom A or B as defined by -m.
# -F: Do not filter data.  See above.
#
# -p: plot.type: one of "depth", "point", or "region"; these require depth, BPC, BPR data files, resp.  "depth" is default.
# -P: Output as PDF file instead of GGP.  This is primarily for convenience and debugging.
# -G fn.ggp: Append graphics to given ggp file, rather than creating new ggp file.
# -B: Annotate for chrom B.  This rotates text and reverses X scale.
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
# * depth: chrom, pos, depth
# * BPC: chromA, posA, chromB, posB, [attribute] 
# * BPR: chromA, posA.start, posA.end, chromB, posB.start, posB.end, [attribute] 
# The attribute column is optional.  See README.file_formats for more details
#
# Read depth may be plotted in several ways:
# * as raw read count 
# * as estimated copy number
# * with log2 scaling
#
# Copy number is estimated as 2 * read depth / RD, where RD is average read depth across the genome,
#   RD = num_reads * read_length / genome_length
# log2 scaling is given as log2(depth/D), where D is median depth of the filtered dataset.
#
# log2 scaling is used when -L flag is set.
# If not log2: if both num.reads and read.length are defined (-u -l) then plot copy number
#   otherwise, plot read depth
#

suppressPackageStartupMessages(library("ggplot2"))
options("width"=300) # change terminal window width
library("reshape2")


set.seed(25)
options(scipen=999)

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
source_relative("BPS_Util.R")
source_relative("DepthUtil.R")

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")

    range.A = parse.range.str(get_val_arg(args, "-A", "all"))  # accessible as range.chr, start=range.pos[1], start=range.pos[2]
    plot.log.depth = get_bool_arg(args, "-L")
    num.reads = as.numeric(get_val_arg(args, "-u", NA))
    read.length = as.numeric(get_val_arg(args, "-l", NA))
    chrom.AB = get_val_arg(args, "-m", "A")
    if (! chrom.AB %in% c("A", "B")) {
        stop("-m argument must be either A or B")
    }

    skip.data.filter = get_bool_arg(args, "-F")
    in.ggp = get_val_arg(args, "-G", NULL)
    plot.type = get_val_arg(args, "-p", "point")
    pdf.out = get_bool_arg(args, "-P")
    is.B = get_bool_arg(args, "-B")

    alpha = get_val_arg(args, "-a", NULL)
    color = get_val_arg(args, "-c", NULL)
    fill = get_val_arg(args, "-f", NULL)
    shape = get_val_arg(args, "-s", NULL)
    linetype = get_val_arg(args, "-t", NULL)
    size = get_val_arg(args, "-z", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)]; args = args[-length(args)]
    data.fn = args[length(args)];      args = args[-length(args)]

    val = list( 
        'range.pos'=range.A$range.pos, 'range.chr'=range.A$range.chr, 'is.B'=is.B, 
        'verbose' = verbose, 'plot.log.depth' = plot.log.depth, 'num.reads' = num.reads,
        'read.length' = read.length, 'skip.data.filter' = skip.data.filter, 'in.ggp' = in.ggp,
        'plot.type' = plot.type, 'pdf.out' = pdf.out, 'alpha' = alpha, 'color' = color, 'fill' = fill,
        'shape' = shape, 'size' = size, 'linetype'=linetype, 'out.ggp' = out.ggp, 'data.fn' = data.fn,
        'chrom.AB' = chrom.AB )
    if (val$verbose) { print(val) }

    return (val)
}

# is.B indicats that this is a panel which will be rotated when assembled in final plot.  
#   * X axis is reversed
make.depth.GGP = function(range.pos = NULL, is.B = FALSE) {
    # define basic properties of Depth GGP object.
    ggp = ggplot() + xlab(NULL) + ylab(NULL) 
    if (!is.null(range.pos)) {
        ggp = ggp + coord_cartesian(xlim=c(range.pos[1], range.pos[2]))
    }
    ggp = ggp + theme_bw()
    ggp = ggp + theme(axis.text=element_text(size=6), axis.text.y=element_text(angle=-90, hjust=0.5))  

    if (is.B) {
        ggp = ggp + scale_x_reverse()  # make x position increase in the downward direction
    }
    # In general, legends will be modified downstream of ggp creation.
    # ggp = ggp + theme(legend.position="none")  
    return(ggp)
}

render.depth = function(ggp, depth, ylabel, alpha=NULL, color=NULL, shape=NULL, size=NULL) {
# if alpha, shape, size, color are defined, their values are passed to corresponding arguments.  
# depth: "chrom", "pos", "depth", "norm.depth"

    # see BreakpointRenderer.R:render.point() for details of how geom_point call constructed
    args = list()   # first collect all static arguments, then add aes and data as arguments
    aes.args = list(x="pos", y="norm.depth")
    if (!is.null(alpha)) args$alpha = alpha
    if (!is.null(shape)) args$shape = shape
    if (!is.null(size)) args$size = size
    if (!is.null(color)) args$color = color
    args$data=depth
    args$mapping = do.call(aes_string, aes.args)  
    # now call geom_point with all the arguments
    ggp = ggp + do.call(geom_point, args)  
    ggp = ggp + ylab(ylabel)
# the simple call
###    ggp = ggp + geom_point(data=depth, aes(x=pos, y=plot.depth), alpha=..., etc.) 
    return(ggp)
}

render.vline = function(ggp, BPC, chrom.AB, alpha=NULL, color=NULL, linetype=NULL, size=NULL) {
# if alpha, size, color, linetype are defined, their values are passed to corresponding arguments.  
# if color is not defined, then color is passed as an aesthetic with data$attribute as value.
# BPC: "chrom.A", "pos.A", "chrom.B", "pos.B", "attribute"
    args = list()   # first collect all static arguments, then add aes and data as arguments
    if (chrom.AB == "A")
        aes.args = list(xintercept="pos.A")
    else
        aes.args = list(xintercept="pos.B")
    if (!is.null(alpha)) args$alpha = alpha
    if (!is.null(linetype)) args$linetype = linetype
    if (!is.null(size)) args$size = size
    if (!is.null(color)) args$color = color
    else { # if not specified, color is an aes with value given by attribute column
        aes.args$color = "attribute"
    }
    args$data=BPC
    args$mapping = do.call(aes_string, aes.args)  
    # now call geom_point with all the arguments
    ggp = ggp + do.call(geom_vline, args)  
# the simple call
# geom_vline(data=rSBP, aes(xintercept=pos), color="#4DAF4A", alpha=0.5) 
    return(ggp)
}

render.CBS = function(ggp, CBS, alpha=NULL, color=NULL, linetype=NULL, size=NULL) {
# if alpha, size, color, linetype are defined, their values are passed to corresponding arguments.  
#  CBS: chrom, start, end, seg.mean, pos, g
    args = list()   
    aes.args = list(x="start", xend="end", y="seg.mean", yend="seg.mean", group="g")
    if (!is.null(alpha)) args$alpha = alpha
    if (!is.null(linetype)) args$linetype = linetype
    if (!is.null(size)) args$size = size
    if (!is.null(color)) args$color = color
    args$data=CBS
    args$mapping = do.call(aes_string, aes.args)  
    ggp = ggp + do.call(geom_segment, args)  
# the simple call
# geom_segment(data=CBS, aes(x=start, xend=end, y=seg.mean, yend=seg.mean, group=g), color="#E41A1C", alpha=1) 
    return(ggp)
}

args = parse_args()

if (is.null(args$in.ggp)) {
    ggp = make.depth.GGP(args$range.pos, args$is.B)
} else {
    ggp = readRDS(args$in.ggp)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
}

if (args$plot.type == "depth" | args$plot.type == "CBS") {
    depth = read.depth(args$data.fn)
    if (!args$skip.data.filter) {
        depth = filter.depth(depth, args$range.chr, args$range.pos)
    }
    data = normalize.depth(depth$depth, args$plot.log.depth, args$num.reads, args$read.length)
    depth$norm.depth = data$norm.depth

    if (args$plot.type == "depth") {
        ggp = render.depth(ggp, depth, data$method, alpha=args$alpha, color=args$color, shape=args$shape, size=args$size)
    } else {  # calculate and overlay CBS
        CBS = getCBS(depth, args$plot.log.depth) 
        ggp = render.CBS(ggp, CBS, alpha=args$alpha, color=args$color, linetype=args$linetype, size=args$size)
    }

} else if (args$plot.type == "point" | args$plot.type == "vline") {
    breakpoints.bpc = read.BPC(args$data.fn)
    if (!args$skip.data.filter) {
        if (args$chrom.AB == "A")
            breakpoints.bpc = filter.BPC(breakpoints.bpc, args$range.chr, args$range.pos, NULL, NULL)
        else
            breakpoints.bpc = filter.BPC(breakpoints.bpc, NULL, NULL, args$range.chr, args$range.pos)
    }
    if (args$plot.type == "point") {
        stop("Unimplemented")
        ggp = render.point(ggp, breakpoints.bpc, args$chrom.AB, alpha=args$alpha, color=args$color, shape=args$shape, size=args$size)
    } else { # (plot.type == "vline")
        ggp = render.vline(ggp, breakpoints.bpc, args$chrom.AB, alpha=args$alpha, color=args$color, linetype=args$linetype, size=args$size)
    }
} else if (args$plot.type == "region" ) {
        stop("Unimplemented")

} else if (args$plot.type == "segment") {
    print("Unimplemented")
} else {
    stop("Unknown plot type", args$plot.type)
}

# now save ggp to file
write.GGP(ggp, args$out.ggp, args$pdf.out)
