# Usage: Rscript HistogramCruncher.R [-v] [-n num.reads] [-l read.length] [-N nbin] [-m hist.max] [-s] [-d]
#       chrom.A.depth.fn chrom.B.depth.fn histogram.ggp
#
# based on breakpointPlotter.R and ReadDepth/CBS/plot_CBS.R

# 6/29/15 - factor of 2 incorporated into copy number calculations

# create a histogram of read depth (possibly normalized to copy number) around integration event for chrom A and B
# These histograms are constructed from read depth files, same format as for read depth plots
# Copy number is obtain by normalizing read depth by (num.reads * read.length) / #bp_in_genome 
# GGP object saved to file given by histogram.ggp
#
# -n num.reads: number of reads in BAM (either total or mapped), typically obtained from .flagstat file
# -l read.length: average read length (from BAM)
# -N nbin: number of bins in histogram.  Default 50
# -m hist.max: max cutoff for copy number (or read depth if calculated) for histogram.
# -v: write verbose output
# -d: draw smooth density plot instead of histogram

suppressPackageStartupMessages(library("ggplot2"))
library("grid")
#library("reshape2")
options("width"=300) # change terminal window width
options(warn=2)

get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments.  Arguments being converted to numbers have NA as default
    verbose = get_bool_arg(args, "-v")
    num.reads = get_val_arg(args, "-n", NA)
    read.length = get_val_arg(args, "-l", NA)
    nbin = get_val_arg(args, "-N", 50)
    hist.max = get_val_arg(args, "-m", NA)
    if (hist.max == "NA") hist.max = NA     # suppress warnings when "NA" is passed to hist.max as an argument
    do.density = get_bool_arg(args, '-d')

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    ggp.fn = args[length(args)];                    args = args[-length(args)]
    virus.depth.all.fn = args[length(args)];        args = args[-length(args)]
    chrom.depth.all.fn = args[length(args)];        args = args[-length(args)]
    chrom.depth.0.fn = args[length(args)];          args = args[-length(args)]

    val = list( 'verbose' = verbose, 'num.reads' = as.numeric(num.reads), 'read.length' = as.numeric(read.length), 'nbin' = as.numeric(nbin),
                'hist.max' = as.numeric(hist.max), 'ggp.fn' = ggp.fn, 'virus.depth.all.fn' = virus.depth.all.fn,
                'chrom.depth.all.fn' = chrom.depth.all.fn, 'chrom.depth.0.fn' = chrom.depth.0.fn, 
                'do.density'=do.density)
    if (val$verbose) { print(val) }
    return (val)
}

estimated.depth = function(num.reads, read.length) {
    num_bp_genome = 3137144693.
    average_depth = num.reads * read.length / num_bp_genome
    return(average_depth)
}

# Read in read depth data and calculate copy number
# read depth is normalized by 2 * (num.reads * read.length) / #bp_in_genome to obtain
#   copy number, with the number of basepairs in genome = 3,137,144,693 
# if num.reads or read.length not specified, copy number is NULL
get.copy.number = function(depth.fn, norm.depth, chrom=NULL, range.start=NA, range.end=NA) {
    read.depth = read.table(depth.fn, header=FALSE, sep="\t", colClasses=c("character","numeric","numeric"), col.names=c("chrom", "pos", "depth"))
    # 14    19000950    21
    # see R/depth_plot.R for more info
    if (!is.na(num.reads) & !is.na(read.length)) {
        read.depth$copy.number = 2. * read.depth$depth / estimated.depth(num.reads, read.length)
    } else {
        cat("Not calculating copy number")
        read.depth$copy.number = NULL
    }
    if (!is.null(chrom)) 
        read.depth = read.depth[read.depth$chrom == chrom,]
    if (!is.na(range.start))
        read.depth = read.depth[read.depth$pos > range.start,]
    if (!is.na(range.end))
        read.depth = read.depth[read.depth$pos <= range.end,]
    return(read.depth)
}

# read in depth data for all files listed in array fn.list
# name.list gives names which identify histogram in legend
# filter.list contains one list for each file, with fields chrom, start, and end
get.histogram.data = function(fn.list, filter.list, name.list, num.reads, read.length) {
    histogram.df = NULL
    # norm.depth is the average read depth per 2 x copy number
    norm.depth = estimated.depth(num.reads, read.length)
    for (i in 1:length(fn.list)) {
        fn=fn.list[i]
        filter = filter.list[[i]]  # filter fields: chrom, start, end
        # pre-calculate average depth and pass that.  Filter fields are optional, skip them.
        d = get.copy.number(fn, norm.depth, filter[[1]], filter[[2]], filter[[3]])
        d$fn = fn
        d$name = name.list[i]
        histogram.df = rbind(histogram.df, d)
    }
    return(histogram.df)
}

# this is specific to virus...
get.fill.scale = function() {
    return(scale_fill_brewer(palette="Set1"))

# if want to set colors explicitly...
#    #Set 1 colors: red, blue, green
#    color.hex = c("#E41A1C", "#377EB8", "#4DAF4A")
#    # want virus red, chrom at breakpoint blue, and chrom all green
#    color.names = c("Virus", "Chrom at breakpoint", "Chrom all")
#    #color.names = c("Chrom at breakpoint", "Chrom all", "Virus")
#    names(color.hex) = color.names
#    return(scale_fill_manual(name="gene.status", values=color.hex))
}

render.histogram = function(histogram.df, nbin, hist.max, do.density) {
    # if copy number has been calculated (because num.reads and read.length specified)
    # then plot that, otherwise plot depth.  We will copy $copy.number to $depth so they're
    # accessed uniformly
    if (!all(is.null(histogram.df$copy.number))) {
        histogram.df$depth = histogram.df$copy.number
        x.label = "Copy Number"
    }  else {
        x.label = "Read Depth"
    }

    # if xmax is a number, limit depth data
    if (!is.na(hist.max)) 
        histogram.df = histogram.df[histogram.df$depth <= hist.max,]

    binwidth = (max(histogram.df$depth) - min(histogram.df$depth)) / nbin

    if (!do.density) {
        p = ggplot(histogram.df, aes(x=depth, fill=name)) + geom_histogram(aes(y=..density..), binwidth=binwidth, alpha=0.6, position="identity") 
    } else {
        p = ggplot(histogram.df, aes(x=depth, fill=name)) + geom_density(alpha=0.6, color=NA)
    }

    p = p + xlab(x.label) + ylab("Frequency") 
    p = p + theme_bw()
    p = p + theme(legend.position=c(0.8,0.8), legend.title=element_text(size=6), legend.text=element_text(size=6), 
            axis.text = element_text(size=6), axis.title=element_text(size=6), legend.key.size=unit(1,"mm"), 
            legend.title=element_blank()) 
    p = p + get.fill.scale()
    return(p)
}

args = parse_args()

# cat(sprintf("Estimated depth %0.2fX\n", estimated.depth(num.reads, read.length)))

# *.list contains data to allow lower-level functions to be agnostic about what they're reading.
# Filter list indicates chromosome and range of depth data to keep (each entry is a triplet, [chrom, start, end] )
#   - Filtering is not used here because anticipate a histogram of entire depth dataset
#if (!args$skip.all.chrom) {  # skip.all.chrom no longer works.
#    fn.list = c(args$chrom.depth.0.fn, args$chrom.depth.all.fn, args$virus.depth.all.fn)
#    #name.list = c(sprintf("Chrom %s at breakpoint", args$chrom), sprintf("Chrom %s all", args$chrom), "Virus")
#    name.list = c("Chrom at breakpoint", "Chrom all", "Virus")
#    filter.list = list( list(args$chrom, args$range.start, args$range.end), list(args$chrom, NA, NA), list(NULL, NA, NA) )
#} else {
fn.list = c(args$chrom.A.depth.fn, args$chrom.B.depth.fn)
name.list = c("Chrom A", "Chrom B")
#filter.list = list( list(args$chrom, args$range.start, args$range.end), list(NULL, NA, NA) )
filter.list = list( list(NULL, NA, NA), list(NULL, NA, NA) )
#}

histogram.df = get.histogram.data(fn.list, filter.list, name.list, args$num.reads, args$read.length)
histogram.ggp = render.histogram(histogram.df, args$nbin, args$hist.max, do.density=args$do.density)

cat(paste("Saving to GPP file", args$ggp.fn, "\n"))
saveRDS(histogram.ggp, file=args$ggp.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/

#ggsave(filename="test.pdf", plot=histogram.ggp, height=6, width=8, useDingbats=FALSE)
# Recall, ggp files can be converted to PDF with ~/bin/ggp2pdf

