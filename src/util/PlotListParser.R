#!/usr/bin/env Rscript

# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Usage: PlotListParser.R [-v] barcode chrom pos PlotList.dat
#
# Given barcode, chrom, and chrom position, return PlotList name which contains this position.
#
# -v: verbose
#

options("width"=180) # useful for debugging

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

# Usage: 
#   args = parse_args()
#   print(args$disease.filter)
parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    # one problem with this approach is that missing mandatory args are hard to catch
    pl.fn = args[length(args)]; args = args[-length(args)]
    pos = as.numeric(args[length(args)]); args = args[-length(args)]
    chrom = args[length(args)]; args = args[-length(args)]
    barcode = args[length(args)]; args = args[-length(args)]

    val = list( 'verbose'=verbose, 'pl.fn' = pl.fn, 'pos' = pos, 'chrom' = chrom, 'barcode' = barcode )
                
    if (val$verbose) { print(val) }

    return (val)
}

args = parse_args()

plot.list = read.table(args$pl.fn, header=TRUE)

#                                 name                      barcode disease virus chrom integration.start integration.end range.start range.end
#1 TCGA-BA-4077-01B-01D-2268-08.chr14A TCGA-BA-4077-01B-01D-2268-08    HNSC HPV16    14          68683065        68742035    68633064  68792035

plot.list = plot.list[plot.list$barcode == args$barcode & plot.list$chrom == args$chrom & 
                      args$pos > plot.list$range.start & args$pos < plot.list$range.end,]

if (nrow(plot.list) == 1) {
    cat(paste(plot.list$name, "\n"))
} else {
    s = sprintf("%s %s:%d", args$barcode, args$chrom, args$pos)

    if (nrow(plot.list) == 0) {
        stop(paste("No lines in", args$pl.fn, "for", s))
    } else {
        stop(paste("Multiple lines in", args$pl.fn, "for", s))
    }
}

