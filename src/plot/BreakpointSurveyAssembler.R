# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript BreakpointSurveyAssembler.R [-v] [-a annotation.A.ggp] [-A annotation.B.ggp] [-H histogram.ggp] [-m] [-N]
#               [-t title] [-h height] [-w width] [-c chrom.A] [-C chrom.B] [-L] [-b] [-d marks.A] [-D marks.B] 
#               breakpoint.ggp depth.A.ggp depth.B.ggp out.pdf
#
# Read in various BreakpointSurveyor GGP objects, assemble them into composite figure, and write to PDF
# in cases where there are no features to annotate (e.g., no genes in range of interest) 
# no annotation.ggp file is passed.  
#
# Input arguments:
# -v: verbose
# -N: Do not align the panels, keeping all ranges and layout as generated (for debugging)
# -L: do not print axis labels
# -b: make big axis text
# -m: format genomic coordinates with commas 
# -d, -D: make alignment marks on chrom A, B resp.  Alignment marks are test lines drawn on various panels to validate that
#         they are aligned correctly.  argument marks.A, marks.B are comma-separated lists genomic position of marks,
#         (e.g., "1000,2000,3000")
# -t: title of plot
# -c, -C: chromosome name (e.g. "11" or "X")
# -h, -w: height, width of PDF output in inches
#
# Layout: 2D breakpoint plot in upper left panel, with Chrom A/B details below/to the left of it, respectively.  Optional
# histogram in lower left.  Chrom details consist of depth and optional annotation.

options("width"=180) # useful for debugging
suppressPackageStartupMessages(library("ggplot2"))
library('grid')
library('gridExtra', quietly=TRUE)
library('gridBase', quietly=TRUE)
library(scales)    # necessary for commas

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

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    annotation.A.fn = get_val_arg(args, "-a", NULL)
    annotation.B.fn = get_val_arg(args, "-A", NULL)
    histogram.fn = get_val_arg(args, "-H", NULL)
    title = get_val_arg(args, "-t", NULL)
    height = as.numeric(get_val_arg(args, "-h", 8))
    width = as.numeric(get_val_arg(args, "-w", 8))
    no.label = get_bool_arg(args, "-L")
    big.font = get_bool_arg(args, "-b")
    chrom.A = get_val_arg(args, "-c", "A")
    chrom.B = get_val_arg(args, "-C", "B")
    commas = get_bool_arg(args, "-m")
    no.align = get_bool_arg(args, "-N")

    marks.A = get_val_arg(args, "-d", NULL)
    marks.B = get_val_arg(args, "-D", NULL)
    if (!is.null(marks.A)) 
        marks.A = unlist(lapply( unlist(strsplit(marks.A, split=",")), as.numeric))
        #marks.A = unlist(lapply( unlist(strsplit(marks.A, split=",")), as.numeric))
    if (!is.null(marks.B)) 
        marks.B = unlist(lapply( unlist(strsplit(marks.B, split=",")), as.numeric))
        #marks.B = unlist(strsplit(marks.B, split=","))

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.fn = args[length(args)];             args = args[-length(args)]
    depth.B.fn = args[length(args)];             args = args[-length(args)]
    depth.A.fn = args[length(args)];             args = args[-length(args)]
    breakpoints.fn = args[length(args)];      args = args[-length(args)]

    val = list('verbose'=verbose, 'breakpoints.fn'= breakpoints.fn, 'histogram.fn'= histogram.fn, 'depth.A.fn'= depth.A.fn,
               'depth.B.fn'= depth.B.fn, 'out.fn'= out.fn, 'annotation.A.fn'=annotation.A.fn, 'no.align'=no.align,
               'annotation.B.fn'=annotation.B.fn, 'title'=title, 'height'=height, 'width'=width, 'marks.A'=marks.A, 'marks.B'=marks.B,
               'no.label'=no.label, 'big.font'=big.font, 'chrom.A'=chrom.A, 'chrom.B'=chrom.B, 'commas'=commas)
    if (val$verbose) { print(val) }

    return (val)
}

# Precise alignment of panels is based on work here: 
    # /Users/mwyczalk/Data/Virus/Virus_2013.9a/RSEM-Exon/RPKM-Scatter/src/RPKM_scatter_plotter.R
# Read more about approach here: http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically

# from http://zevross.com/blog/2014/11/20/under-the-hood-of-ggplot2-graphics-in-r/
    #The ggplot_build function outputs a list of data frames (one for each layer)
    #and a panel object with information about axis limits among other things. The
    #ggplot_gtable function, which takes the ggplot_build object as input, builds
    #all grid graphical objects (known as “grobs”) necessary for displaying the
    #plot. You can manipulate the output from ggplot_build.
# The general idea is to set the ranges of .gb objects first, then the widths/heights of .gt objects
# details of .gb and .gt objects can be obtained with print(str(foo))
# .ggp and .gt objects can be visualized with grid.draw(foo)

# set range and ticks in X direction.  Return target.gb
align_X_range = function(target.gb, template.gb) {
    target.gb$panel$ranges[[1]]$x.range = template.gb$panel$ranges[[1]]$x.range
    target.gb$panel$ranges[[1]]$x.labels = template.gb$panel$ranges[[1]]$x.labels
    target.gb$panel$ranges[[1]]$x.major = template.gb$panel$ranges[[1]]$x.major
    target.gb$panel$ranges[[1]]$x.minor = template.gb$panel$ranges[[1]]$x.minor
    target.gb$panel$ranges[[1]]$x.major_source = template.gb$panel$ranges[[1]]$x.major_source
    target.gb$panel$ranges[[1]]$x.minor_source = template.gb$panel$ranges[[1]]$x.minor_source
    return(target.gb)
}

# set range and ticks in Y direction.  Return target.gb
align_Y_range = function(target.gb, template.gb) {
    target.gb$panel$ranges[[1]]$y.range = template.gb$panel$ranges[[1]]$y.range
    target.gb$panel$ranges[[1]]$y.labels = template.gb$panel$ranges[[1]]$y.labels
    target.gb$panel$ranges[[1]]$y.major = template.gb$panel$ranges[[1]]$y.major
    target.gb$panel$ranges[[1]]$y.minor = template.gb$panel$ranges[[1]]$y.minor
    target.gb$panel$ranges[[1]]$y.major_source = template.gb$panel$ranges[[1]]$y.major_source
    target.gb$panel$ranges[[1]]$y.minor_source = template.gb$panel$ranges[[1]]$y.minor_source
    return(target.gb)
}

# align the range and tick marks of target object to match template object.  Dim is one of X, Y, or XY.
# returns target.gb, a ggplot_build object.  Template is unchanged.
# if no.align = true, convert target.ggp to target.gb, but don't modify range
align_range = function(target.ggp, template.ggp, dim, no.align) {
    # The ggplot_build function outputs a list of data frames (one for each layer)
    # and a panel object with information about axis limits among other things.
    # You can manipulate the output from ggplot_build.  $panel is the list element we're interested in.
    target.gb = ggplot_build(target.ggp)
    template.gb = ggplot_build(template.ggp)

    if (!no.align) {
        # Setting the range will subset and scale the data. Set tick marks explitly.  
        if (dim == "X") {
            target.gb = align_X_range(target.gb, template.gb)
        } else if (dim == "Y") {
            target.gb = align_Y_range(target.gb, template.gb)
        } else if (dim == "XY") {
            target.gb = align_X_range(target.gb, template.gb)
            target.gb = align_Y_range(target.gb, template.gb)
        } else {
            stop("Unknown dim", dim)
        }
    }
    return(target.gb)
}

# convert x .gt objects (gtable) object.
# x can be a ggp (ggplot), gb (ggplot_build), or gt object
# return x.gt object
convert.gt = function(x) {
    if (inherits(x, "ggplot")) { # is ggp
        x.gt = ggplot_gtable(ggplot_build(x))
    } else if (inherits(x, "gtable")) {  # is gt
        x.gt = x
    } else if (is.list(x)) {  # assuming .gb
        x.gt = ggplot_gtable(x)
    } else stop("Unknown type")
    return(x.gt)
}

# set the width and/or height of target to match template.
# target and template can be .ggp, .gb, or .gt objects
# dim is one of X, Y, XY to make target equal the width, height, or both of template
# return target.gt
align_layout = function(target, template, dim, no.align) {
    # first, convert target and template to .gt objects, whatever they may currently be
    target.gt = convert.gt(target)
    template.gt = convert.gt(template)
    
    if (!no.align) {
        # (The ggplot_gtable function, which takes the ggplot_build object as input,
        # builds all grid graphical objects (known as “grobs”) necessary for displaying
        # the plot.)
        # now match the widths and/or the heights
        if (dim == "X") {
            target.gt$widths = template.gt$widths
        } else if (dim == "Y") {
            target.gt$heights = template.gt$heights
        } else if (dim == "XY") {
            target.gt$widths = template.gt$widths
            target.gt$heights = template.gt$heights
        } else {
            stop("Unknown dim", dim)
        }
    }
    return(target.gt)
}

# This is useful for debugging of data range and tick mark position.  
add.alignment.marks = function(ggp, chrom.marks.v=NULL, chrom.marks.h=NULL) {
    if (!is.null(chrom.marks.v))
        ggp = ggp + geom_vline(xintercept=chrom.marks.v)
    if (!is.null(chrom.marks.h))
        ggp = ggp + geom_hline(yintercept=chrom.marks.h)
    return(ggp)
}

# return target.gt which has both range and layout matching that of template along given axis ("X" or "Y")
# if no.align is TRUE, skip alignment but do all other conversions
align.panels = function(target.ggp, template.ggp, axis, no.align=FALSE) {
    if (! axis %in% c("X", "Y")) stop("Unknown axis", axis)
    # first set the data range and labels
    target.gb = align_range(target.ggp, template.ggp, axis, no.align)
    # next set the panel widths/heights
    target.gt = align_layout(target.gb, template.ggp, axis, no.align)
    return(target.gt)
}

# This figure includes only breakpoint plot and read depth.
# if no.align = TRUE, assemble panels but do not align them (change their range or layout)
assemble_plot = function(breakpoint.ggp, depth.A.ggp, depth.B.ggp, histogram.ggp, annotation.A.ggp, annotation.B.ggp, title.ggp, no.label=FALSE, no.align=FALSE) {

    # get rid of y-axis titles on depth
    depth.A.ggp = depth.A.ggp + theme(axis.title.y = element_blank())
    depth.B.ggp = depth.B.ggp + theme(axis.title.x = element_blank())

    depth.A.gt = align.panels(depth.A.ggp, breakpoint.ggp, "X", no.align)
    depth.B.gt = align.panels(depth.B.ggp, breakpoint.ggp, "Y", no.align)

    # make composite depth + annotation panels as necessary
    if (!is.null(annotation.A.ggp)) {
        annotation.A.gt = align.panels(annotation.A.ggp, breakpoint.ggp, "X", no.align)
        panel.A.gt = arrangeGrob(annotation.A.gt, depth.A.gt, heights=c(0.5,0.5), ncol=1, nrow=2)
    } else {
        panel.A.gt = depth.A.gt
    }

    if (!is.null(annotation.B.ggp)) {
        annotation.B.ggp = annotation.B.ggp + coord_flip()
        annotation.B.gt = align.panels(annotation.B.ggp, breakpoint.ggp, "Y", no.align)
        panel.B.gt = arrangeGrob(depth.B.gt, annotation.B.gt, widths=c(0.5,0.5), ncol=2, nrow=1)
    } else {
        panel.B.gt = depth.B.gt
    }

    if (is.null(histogram.ggp)) {
        histogram.ggp = grid.rect(gp=gpar(col="white"))
        #histogram.ggp = rectGrob(gp = gpar(col = "white"))  # this may not be necessary
    } 
    # not sure this is necessary
#    else {
    # make sure histogram dimensions match
##        histogram.gt = ggplot_gtable(ggplot_build(histogram.ggp))
##        histogram.gt$heights = depth.A.gt$heights
##        histogram.gt$widths = depth.B.gt$widths
#        histogram.gt = align_layout(histogram.ggp, depth.A.gt, "Y")
#        histogram.gt = align_layout(histogram.gt, depth.B.gt, "X")
#    }

    main.grob = arrangeGrob(panel.B.gt, breakpoint.ggp, histogram.ggp, panel.A.gt, widths=c(0.3,0.7), heights=c(0.7,0.3), ncol=2, nrow=2)

    if (!is.null(title.ggp)) {
        annotated.grob = grid.arrange(title.ggp, main.grob, heights = c(0.025, 0.975), ncol=1, nrow=2, newpage=FALSE)
    } else {
        grid.draw(main.grob)
    }

#    if (!no.label) { 
#        # hard positioning of labels.  Ick.
#        grid.text("Chromosome Copy Number", x = unit(0.160, "npc"), y = unit(0.285, "npc"), gp=gpar(fontsize=10))
#        grid.text("Virus Copy Number", x = unit(0.295, "npc"), y = unit(0.150, "npc"), gp=gpar(fontsize=10), rot=-90)
#    }
}


args = parse_args()

breakpoint.ggp = readRDS(args$breakpoints.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/

depth.A.ggp = readRDS(args$depth.A.fn)
depth.B.ggp = readRDS(args$depth.B.fn)

annotation.A.ggp = NULL
annotation.B.ggp = NULL
histogram.ggp = NULL
title.ggp = NULL

if (!is.null(args$annotation.A.fn))
    annotation.A.ggp = readRDS(args$annotation.A.fn)

if (!is.null(args$annotation.B.fn))
    annotation.B.ggp = readRDS(args$annotation.B.fn)
if (!is.null(args$histogram.fn)) 
    histogram.ggp = readRDS(args$histogram.fn)

title_text = args$title
if (!is.null(title_text))  {
    title.ggp = textGrob(title_text)
} else {
    title.ggp = NULL
}

gray.ticks = theme(axis.ticks = element_line(color="gray50"))
no.margin = theme(plot.margin=unit(c(0,0,0,0),"in"))

breakpoint.ggp = breakpoint.ggp + theme(plot.margin=unit(c(0.125,0.125,0,0), "in")) # top, right, bottom, and left margins
depth.A.ggp = depth.A.ggp + xlab(sprintf("Chr %s Position", args$chrom.A)) + no.margin
depth.B.ggp = depth.B.ggp + theme( axis.text.y=element_text(angle=-90, hjust=0.5), axis.title.y = element_text(angle=-90)) + xlab(sprintf("Chr %s Position", args$chrom.B)) + coord_flip() + no.margin
annotation.A.ggp = annotation.A.ggp + no.margin
annotation.B.ggp = annotation.B.ggp + no.margin

if (args$no.label & !is.null(histogram.ggp)) {
    histogram.ggp = histogram.ggp + theme(axis.title = element_blank(), legend.position="none")
}

# commas on genomic positions.  Implmement this somewhere
# from ~/Data/TCGA_SARC/ICGC/BreakpointSurveyor/E_Breakpoint/src/BreakpointRenderer.R
#    if (commas) {
#        p = p + scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma)
#    }

if (!is.null(args$marks.A)) {
    depth.A.ggp = add.alignment.marks(depth.A.ggp, args$marks.A)
    breakpoint.ggp = add.alignment.marks(breakpoint.ggp, args$marks.A, NULL)
    annotation.A.ggp = add.alignment.marks(annotation.A.ggp, args$marks.A)
}
if (!is.null(args$marks.B)) {
    depth.B.ggp = add.alignment.marks(depth.B.ggp, args$marks.B)
    breakpoint.ggp = add.alignment.marks(breakpoint.ggp, NULL, args$marks.B)
    annotation.B.ggp = add.alignment.marks(annotation.B.ggp, args$marks.B)
}

#    chrom.A.marks = c(33080000, 33140000)
#    chrom.B.marks = c(120825000, 120900000)
#    depth.A.ggp = add.alignment.marks(depth.A.ggp, chrom.A.marks)
#    depth.B.ggp = add.alignment.marks(depth.B.ggp, chrom.B.marks)
#    breakpoint.ggp = add.alignment.marks(breakpoint.ggp, chrom.A.marks, chrom.B.marks)
#    annotation.A.ggp = add.alignment.marks(annotation.A.ggp, chrom.A.marks)
#    annotation.B.ggp = add.alignment.marks(annotation.B.ggp, chrom.B.marks)


# we're expanding meaning of big.font to do many customizaitons for -special figures for manuscript
# This should be done in an external script
if (args$big.font) {
    target.size = 12
    no.grid = theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
    resize.text = theme(axis.text = element_text(size=target.size, color="gray50"))
    if (!is.null(histogram.ggp)) histogram.ggp = histogram.ggp + no.grid + resize.text + gray.ticks
    depth.A.ggp = depth.A.ggp + no.grid + resize.text + gray.ticks
    depth.B.ggp = depth.B.ggp + no.grid + resize.text + gray.ticks
    breakpoint.ggp = breakpoint.ggp + no.grid + theme(axis.text = element_text(size=6, color="gray50")) + gray.ticks

    if (!is.null(histogram.ggp)) histogram.ggp = histogram.ggp + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

cat(sprintf("Saving to %s\n", args$out.fn))
pdf(file=args$out.fn, width=args$width, height=args$height, useDingbats=FALSE)
assemble_plot(breakpoint.ggp, depth.A.ggp, depth.B.ggp, histogram.ggp, annotation.A.ggp, annotation.B.ggp, title.ggp, args$no.label, args$no.align)

