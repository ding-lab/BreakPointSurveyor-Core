# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Usage: Rscript DepthDrawer.R [-v] [-P] [-M range] [-N range] [-F] [-G fn.ggp] [-p plot.type] [-l]
#                [-u num.reads] [-n read.length] [-m chrom] [-L] [-B] [-b] [-C] [-y y.mid] [-j y.jitter]
#                [-a alpha] [-c color] [-f fill] [-s shape] [-z size] [-t linetype] data.fn out.ggp
#
#   Plot read depth (or related quantites) over a genomic region and add annotation to this plot.
#
#   These values of plot.type (which define the data types to be drawn or added) are supported:
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
# -n: read length.  Necessary for normalizing to copy number
#
# -m chrom: if plotting point or region from BPC or BPR data, resp., define whether using chrom A or B.  "A" is default.
# -M: Define range for this chromosome, specified as "C" or "C:r-s", where C is chromosome name
#         and r,s are genomic start, stop positions.  Range is used for two distinct purposes:
#           1) Filter data before plotting
#              Filtering data may be prevented with -F.  Applies to chrom A or B as defined by -m.
#           2) Specify plot region
#              Plot region is specified only if -G is not specified (i.e., creating new GGP)
# -N Define range for plotting opposite chromosome.  Defines filtering of BPC/BPR data before plotting.
#    See below for details.
# -C: Use value of "other" (chrom N) chromosome as color for BPC/BPR-based annotation.  Attributes in those files or -c 
#     will override this.
# -F: Do not filter data for chrom M.  See above.
# -l: Swap A, B columns when reading in BPC or BPR, so chrom A > chrom B by string comparison
#
# -p: plot.type: one of "depth", "point", or "region"; these require depth, BPC, BPR data files, resp.  "depth" is default.
# -P: Output as PDF file instead of GGP.  This is primarily for convenience and debugging.
# -G fn.ggp: Append graphics to given ggp file, rather than creating new ggp file.
# -B: Annotate for panel B in assembled figure  This rotates text and reverses X scale.
# -b: format genomic coordinates without commas 
#
# If plot.type is "point", vertical position is random within y.mid +/- y.jitter.  If not defined, both values
# taken from range of plot
# -y y.mid: define y.mid
# -j y.jitter: define y.jitter
#
# The following constant attributes can be defined:
#   -a alpha
#   -c color  [e.g., "-c gray10" or "-c #FF0000"]
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
# If not log2: if both num.reads and read.length are defined (-u -n) then plot copy number
#   otherwise, plot read depth
#
# Specifying chrom and range for BPC/BPR files (TODO: simplify and integrate details below)
# When plotting read depth files, data is specific to that chrom, and -M defines the range of interest.
# However, BPC (and BPR) files have chrom A and B, and we must specify which to plot.
#   -m specifies whether chrom associated with read depth is chrom A or B.
# Then, chrom.M is the one whose depth we are plotting, and chrom.N is the mate chrom
#   (i.e., -m B implies that chrom.M is B and chrom.N is A)
# When plotting BPC/BPR data as annotation to read depth, 
#   -M specifies the chrom name and range of chrom.M (e.g., chr1:123-456)
#   -N specifies the filter parameter for the mate chrom.  Specifying both chrom name and range
#      will include only those breakpoints plotted in Breakpoints Coordinate panel.  Alternatively,
#      excluding chrom or range specification will include chrom.B breakpoints from all chrom and all
#      genomic positions, resp.  e.g., ("-m A -M chr1:100-200 -N 500-600" will annotate read depth plot with 
#      marks indicating the position along chrom.A ("chr1") of all breakpoints (position 100-200) for which breakpoints along chrom B 
#      occur at positions 500-600 on any chrom.
#      Note that default value of -N will not exclude any chrom.N breakpoints.

suppressPackageStartupMessages(library("ggplot2"))
options("width"=300) # change terminal window width
library("reshape2")
library(scales)    # necessary for commas


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
source_relative("../util/BPS_Util.R")
source_relative("BPS_PlotUtil.R")
source_relative("DepthUtil.R")

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")

    range.M = parse.range.str(get_val_arg(args, "-M", NULL))  
    range.N = parse.range.str(get_val_arg(args, "-N", NULL))  
    plot.log.depth = get_bool_arg(args, "-L")
    num.reads = as.numeric(get_val_arg(args, "-u", NA))
    read.length = as.numeric(get_val_arg(args, "-n", NA))
    chrom.AB = get_val_arg(args, "-m", "A")
    if (! chrom.AB %in% c("A", "B")) {
        stop("-m argument must be either A or B")
    }
    if (chrom.AB == "A") {
        range.A = range.M
        range.B = range.N
    } else {
        range.B = range.M
        range.A = range.N
    }

    skip.data.filter = get_bool_arg(args, "-F")
    in.ggp = get_val_arg(args, "-G", NULL)
    plot.type = get_val_arg(args, "-p", "point")
    pdf.out = get_bool_arg(args, "-P")
    panel.B = get_bool_arg(args, "-B")
    no.commas = get_bool_arg(args, "-b")
    color.by.chrom.N = get_bool_arg(args, "-C")
    flip.ab = get_bool_arg(args, "-l")

    alpha = as.numeric(get_val_arg(args, "-a", NA))
    color = get_val_arg(args, "-c", NULL)
    fill = get_val_arg(args, "-f", NULL)
    shape = as.numeric(get_val_arg(args, "-s", NA))
    linetype = get_val_arg(args, "-t", NULL)
    size = as.numeric(get_val_arg(args, "-z", NA))

    y.mid = as.numeric(get_val_arg(args, "-y", NA))
    y.jitter = as.numeric(get_val_arg(args, "-j", NA))

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)]; args = args[-length(args)]
    data.fn = args[length(args)];      args = args[-length(args)]

    val = list( 
        'range.A.pos'=range.A$range.pos, 'range.A.chr'=range.A$range.chr, 
        'range.B.pos'=range.B$range.pos, 'range.B.chr'=range.B$range.chr, 
        'range.M.pos'=range.M$range.pos, 'range.M.chr'=range.M$range.chr, 
        'range.N.pos'=range.N$range.pos, 'range.N.chr'=range.N$range.chr, 
        'panel.B'=panel.B, 'y.mid'=y.mid, 'y.jitter'=y.jitter,
        'verbose' = verbose, 'plot.log.depth' = plot.log.depth, 'num.reads' = num.reads,
        'read.length' = read.length, 'skip.data.filter' = skip.data.filter, 'in.ggp' = in.ggp,
        'plot.type' = plot.type, 'pdf.out' = pdf.out, 'alpha' = alpha, 'color' = color, 'fill' = fill,
        'shape' = shape, 'size' = size, 'linetype'=linetype, 'out.ggp' = out.ggp, 'data.fn' = data.fn,
        'chrom.AB' = chrom.AB, 'no.commas'=no.commas, 'color.by.chrom.N'=color.by.chrom.N, 'flip.ab'=flip.ab )
    if (val$verbose) { print(val) }

    return (val)
}

# panel.B indicates that this is a panel which will be rotated when assembled in final plot, and we reverse X axis
# (making genomic position increase in downward direction in final plot)
make.depth.GGP = function(range.pos = NULL, panel.B = FALSE, no.commas=FALSE) {
    # define basic properties of Depth GGP object.
    ggp = ggplot() + xlab(NULL) + ylab(NULL) 
    if (!is.null(range.pos)) {
        ggp = ggp + coord_cartesian(xlim=c(range.pos[1], range.pos[2]))
    }
    ggp = ggp + theme_bw()
    ggp = ggp + theme(axis.text=element_text(size=6), axis.text.y=element_text(angle=-90, hjust=0.5))  

    trans = if (panel.B) reverse_trans() else "identity"
    labels = if (no.commas) waiver() else comma
    ggp = ggp + scale_x_continuous(labels=labels, trans=trans)  
        
    return(ggp)
}

render.depth = function(ggp, depth, ylabel, alpha=NA, color=NULL, shape=NA, size=NA) {
# if alpha, shape, size, color are defined, their values are passed to corresponding arguments.  
# depth: "chrom", "pos", "depth", "norm.depth"

    # see BreakpointDrawer.R:render.point() for details of how geom_point call constructed
    args = list()   # first collect all static arguments, then add aes and data as arguments
    aes.args = list(x="pos", y="norm.depth")
    if (!is.na(alpha)) args$alpha = alpha
    if (!is.na(shape)) args$shape = shape
    if (!is.na(size)) args$size = size
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

# draw points at X position given by BPC$pos and y position "jittered", placed randomly at y.mid +/- y.jitter
render.point = function(ggp, BPC, chrom.AB, y.mid, y.jitter, alpha=NA, color=NULL, fill=NULL, shape=NA, size=NA) {
# BPC: "chrom.A", "pos.A", "chrom.B", "pos.B", "attribute"
    # see BreakpointDrawer.R:render.point() for details of how geom_point call constructed
    args = list()   # first collect all static arguments, then add aes and data as arguments

    if (chrom.AB == "A")
        aes.args = list(x="pos.A", y=y.mid)
    else
        aes.args = list(x="pos.B", y=y.mid)
    args$position = position_jitter(height=y.jitter, width=0)

    if (!is.na(alpha)) args$alpha = alpha
    if (!is.na(shape)) args$shape = shape
    if (!is.na(size)) args$size = size
    if (!is.null(color)) args$color = color
    if (!is.null(fill)) args$fill = fill
    else { # if not specified, color is an aes with value given by attribute column
        aes.args$color = "attribute"
    }
    args$data=BPC
    args$mapping = do.call(aes_string, aes.args)  
    # now call geom_point with all the arguments
    ggp = ggp + do.call(geom_point, args)  
# the simple call
###    ggp = ggp + geom_point(data=depth, aes(x=pos, y=plot.depth), alpha=..., etc.) 
    return(ggp)
}

render.vline = function(ggp, BPC, chrom.AB, alpha=NA, color=NULL, linetype=NULL, size=NA) {
# if alpha, size, color, linetype are defined, their values are passed to corresponding arguments.  
# if color is not defined, then color is passed as an aesthetic with data$attribute as value.
# BPC: "chrom.A", "pos.A", "chrom.B", "pos.B", "attribute"
    args = list()   # first collect all static arguments, then add aes and data as arguments
    if (chrom.AB == "A")
        aes.args = list(xintercept="pos.A")
    else
        aes.args = list(xintercept="pos.B")
    if (!is.na(alpha)) args$alpha = alpha
    if (!is.null(linetype)) args$linetype = linetype
    if (!is.na(size)) args$size = size
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

# Render region based on BPR data
# names(BPR) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "attribute")
render.region = function(ggp, BPR, chrom.AB, alpha=NA, color=NULL, fill=NULL) {
    # if alpha, fill, or color are defined, their values are passed to corresponding arguments.  
    # if fill is not defined, then fill is passed as an aesthetic with data$attribute as value.
    if (is.null(BPR)) return (ggp)

    args = list()   # first collect all static arguments, then add aes and data as arguments
    if (chrom.AB == "A") {
        aes.args = list(xmin="pos.A.start", xmax="pos.A.end", ymin=-Inf, ymax=Inf)
    } else {
        aes.args = list(xmin="pos.B.start", xmax="pos.B.end", ymin=-Inf, ymax=Inf)
    }

    if (!is.na(alpha)) args$alpha = alpha
    if (!is.null(color)) args$color = color
    if (!is.null(fill)) args$fill = fill
    else { # if not specified, fill is an aes with value given by attribute column
        aes.args$fill = "attribute"
    }
    args$data=BPR
    args$mapping = do.call(aes_string, aes.args)  
    ggp = ggp + do.call(geom_rect, args)  

    # If there was no attribute column specified in BPR file (all attributes are "default") then
    # get rid of color legend.
#    attrib.uniq = unique(BPR$attribute)
#    if (length(attrib.uniq) == 1 & attrib.uniq[1] == "default") {
#        ggp = ggp + scale_fill_discrete(guide=FALSE)  
#    }

    return(ggp)
}

render.CBS = function(ggp, CBS, alpha=NA, color=NULL, linetype=NULL, size=NA) {
# if alpha, size, color, linetype are defined, their values are passed to corresponding arguments.  
#  CBS: chrom, start, end, seg.mean, pos, g
    args = list()   
    aes.args = list(x="start", xend="end", y="seg.mean", yend="seg.mean", group="g")
    if (!is.na(alpha)) args$alpha = alpha
    if (!is.null(linetype)) args$linetype = linetype
    if (!is.na(size)) args$size = size
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
    ggp = make.depth.GGP(args$range.pos, args$panel.B, args$no.commas)
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
    breakpoints = read.BPC(args$data.fn, args$flip.ab)
    if (!args$skip.data.filter) {
        breakpoints = filter.BPC(breakpoints, args$range.A.chr, args$range.A.pos, args$range.B.chr, args$range.B.pos)
    }

    if (nrow(breakpoints) > 0) {
        if (args$color.by.chrom.N & all(breakpoints$attribute == "default")) {
            if (args$chrom.AB == "A")
                breakpoints$attribute = breakpoints$chrom.B
            else
                breakpoints$attribute = breakpoints$chrom.A
        }

        if (args$plot.type == "point") {  
            # Points have a random Y position y.pos +/- y.jitter
            # obtain y.pos and jitter range from existing range of plot
            # This can be optionally specified at command line
            # See BreakpointSurveyAssembler.R:align_X_range() for details about accessing range.
            #   requires ggplot2 >= 2.2.0

            y.range = ggplot_build(ggp)$layout$panel_ranges[[1]]$y.range
            y.jitter = (y.range[2] - y.range[1]) / 2
            y.mid = y.range[1] + y.jitter

            y.mid = if (is.na(args$y.mid)) y.mid else args$y.mid
            y.jitter = if (is.na(args$y.jitter)) y.jitter else args$y.jitter

            ggp = render.point(ggp, breakpoints, args$chrom.AB, y.mid, y.jitter, 
                alpha=args$alpha, color=args$color, shape=args$shape, size=args$size)
        } else { # (plot.type == "vline")
            ggp = render.vline(ggp, breakpoints, args$chrom.AB, alpha=args$alpha, color=args$color, linetype=args$linetype, size=args$size)
        }
    }
} else if (args$plot.type == "region" ) {
    breakpoints = read.BPR(args$data.fn, args$flip.ab)
    if (!args$skip.data.filter) {
        breakpoints = filter.BPR(breakpoints, args$range.A.chr, args$range.A.pos, args$range.B.chr, args$range.B.pos)
    }

    if (args$color.by.chrom.N & all(breakpoints$attribute == "default")) {
        if (args$chrom.AB == "A")
            breakpoints$attribute = breakpoints$chrom.B
        else
            breakpoints$attribute = breakpoints$chrom.A
    }

    ggp = render.region(ggp, breakpoints, args$chrom.AB, alpha=args$alpha, color=args$color, fill=args$fill)

} else if (args$plot.type == "segment") {
    print("Unimplemented")
} else {
    stop("Unknown plot type", args$plot.type)
}

# now save ggp to file
write.GGP(ggp, args$out.ggp, args$pdf.out)
