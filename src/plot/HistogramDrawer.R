# Usage: Rscript HistogramDrawer.R [-v] [-u num.reads] [-n read.length] [-N nbin] [-m hist.max] [-d] [-P]
#       [-e names] depth.A.fn depth.B.fn out.ggp
#
# create a histogram of read depth (or estimated copy number) for chrom A and B
# These histograms are constructed from read depth data
# Copy number is obtained by normalizing read depth by (num.reads * read.length) / 2 * #bp_in_genome 
# GGP object saved to file given by out.ggp
#
# -u num.reads: number of reads in BAM (either total or mapped), typically obtained from .flagstat file
#               Necessary for normalizing to copy number
# -n read.length: average read length (from BAM)
# -N nbin: number of bins in histogram.  Default 50
# -m hist.max: max cutoff for read depth / copy number for histogram.
# -e names: comma-separated list of labels for legend
# -v: write verbose output
# -d: draw smooth density plot instead of histogram
# -P: Output as PDF file instead of GGP.  This is primarily for convenience and debugging.

suppressPackageStartupMessages(library("ggplot2"))
library("grid")
#library("reshape2")
options("width"=300) # change terminal window width
options(warn=2)


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

    # optional arguments.  Arguments being converted to numbers have NA as default
    verbose = get_bool_arg(args, "-v")
    num.reads = get_val_arg(args, "-u", NA)
    read.length = get_val_arg(args, "-n", NA)
    legend.names = unlist(strsplit(get_val_arg(args, "-e", "Chrom A,Chrom B"), split=","))
    nbin = get_val_arg(args, "-N", 50)
    hist.max = get_val_arg(args, "-m", NA)
    #if (hist.max == "NA") hist.max = NA     # suppress warnings when "NA" is passed to hist.max as an argument
    do.density = get_bool_arg(args, '-d')
    pdf.out = get_bool_arg(args, '-P')

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)];           args = args[-length(args)]
    depth.B.fn = args[length(args)];        args = args[-length(args)]
    depth.A.fn = args[length(args)];        args = args[-length(args)]

    val = list( 'verbose' = verbose, 'num.reads' = as.numeric(num.reads), 'read.length' = as.numeric(read.length), 'nbin' = as.numeric(nbin),
                'hist.max' = as.numeric(hist.max), 'out.ggp' = out.ggp, 'depth.A.fn' = depth.A.fn, 'pdf.out'=pdf.out,
                'depth.B.fn' = depth.B.fn, 'do.density'=do.density, 'legend.names'=legend.names)
    if (val$verbose) { print(val) }
    return (val)
}

# Read in depth data for a set of files, optionally filter and normalize, and combine into
# single dataset.
# data.list contains information about all datasets to be read: (DS1, DS2, ...)
#   each DS is a list (fn, name, filter)
#     fn is filename of depth dataset 
#     name is human-readable identifier of dataset
#     filter is list of (range.chr, range.pos), where range.pos is list (start, end)
#       filter may be NULL to skip filtering
# return list (histogram.df, method) where method is a text string indicating "copy number" or "read depth"
get.histogram.data = function(data.list, num.reads=NULL, read.length=NULL) {
    histogram.df = NULL
    for (i in 1:length(data.list)) {
        ds=data.list[[i]]
        depth = read.depth(ds$fn)
        if (!is.null(ds$filter)) {
            depth = filter.depth(depth, fn$filter[1], fn$filter[2])
        }
        data = normalize.depth(depth$depth, FALSE, num.reads, read.length)
        depth$norm.depth = data$norm.depth
        depth$name = ds$name
        histogram.df = rbind(histogram.df, depth)
    }
    # relying on fact that method is same for all DS
    return(list("histogram.df"=histogram.df, "method"=data$method))
}

get.fill.scale = function() {
    return(scale_fill_brewer(palette="Set1"))
}

render.histogram = function(histogram.df, x.label, nbin, hist.max, do.density) {
    # if xmax is a number, limit depth data
    if (!is.na(hist.max)) 
        histogram.df = histogram.df[histogram.df$norm.depth <= hist.max,]

    binwidth = (max(histogram.df$norm.depth) - min(histogram.df$norm.depth)) / nbin

    if (!do.density) {
        ggp = ggplot(histogram.df, aes(x=norm.depth, fill=name)) + geom_histogram(aes(y=..density..), binwidth=binwidth, alpha=0.6, position="identity") 
    } else {
        ggp = ggplot(histogram.df, aes(x=norm.depth, fill=name)) + geom_density(alpha=0.6, color=NA)
    }

    ggp = ggp + xlab(x.label) + ylab("Frequency") 
    ggp = ggp + theme_bw()
    ggp = ggp + theme(legend.position=c(0.8,0.8), legend.title=element_text(size=6), legend.text=element_text(size=6), 
            axis.text = element_text(size=6), axis.title=element_text(size=6), legend.key.size=unit(1,"mm"))
    ggp = ggp + theme(legend.title=element_blank()) 
    ggp = ggp + get.fill.scale()
    return(ggp)
}

args = parse_args()

# data.list defines the filename, label, and filter (if any) of each dataset in histogram.
# currently we only have two datasets (depth A and depth B) and these are not being filtered.
data.list = list(  list(fn=args$depth.A.fn, name=args$legend.names[1], filter=NULL), list(fn=args$depth.B.fn, name=args$legend.names[2], filter=NULL))

data = get.histogram.data(data.list, args$num.reads, args$read.length)
ggp = render.histogram(data$histogram.df, data$method, args$nbin, args$hist.max, do.density=args$do.density)

write.GGP(ggp, args$out.ggp, args$pdf.out)
