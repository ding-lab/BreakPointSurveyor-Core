#!/usr/bin/env Rscript

# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Usage: ZoomGGP.R [-v] [-x X] [-y Y] [-r r.x] [-R r.y] [-h height] [-w width] [-F] [-s x,y] data.ggp out.pdf
#
# Change plot limits of a GGP file and save as PDF
#
# -v: verbose
# -h: height in inches of output PDF (default 6)
# -w: width in inches of output PDF (default 8)
# -x: center X-axis of plot at ths value
# -y: center Y-axis of plot at ths value
# -r: set extents ('radius') in X direction
# -Y: set extents ('radius') in Y direction
# -F: do not swap X,Y coords.  (note that discordant plots do coord_flip() which does swap X,Y)
# -s x,y: draw a "special" mark at coordinates given by x,y (comma-delimited pair of numbers)
#

suppressPackageStartupMessages(library("ggplot2"))

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
    height = as.numeric(get_val_arg(args, "-h", "6"))
    width = as.numeric(get_val_arg(args, "-w", "8"))
    x.mid = as.numeric(get_val_arg(args, "-x", NA))
    y.mid = as.numeric(get_val_arg(args, "-y", NA))
    x.rad = as.numeric(get_val_arg(args, "-r", NA))
    y.rad = as.numeric(get_val_arg(args, "-R", NA))
    no.flip = get_bool_arg(args, "-F")
    special.xy = get_val_arg(args, "-s", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    # one problem with this approach is that missing mandatory args are hard to catch
    out.fn = args[length(args)]; args = args[-length(args)]
    ggp.fn = args[length(args)]; args = args[-length(args)]

    # parse out x,y values of special.xy, if passed
    if (!is.null(special.xy)) {
        xy=strsplit(special.xy, ',', fixed=TRUE)[[1]]
        special.x = as.numeric(xy[1])
        special.y = as.numeric(xy[2])
    } else {
        special.x = NA
        special.y = NA
    }

    val = list( 'verbose'=verbose, 'height'=as.numeric(height), 'width'=as.numeric(width),
                'x.mid'=x.mid, 'y.mid'=y.mid, 'x.rad'=x.rad, 'y.rad'=y.rad, 'no.flip'=no.flip,
                'out.fn'=out.fn, 'ggp.fn'=ggp.fn, 'special.x'=special.x, 'special.y'=special.y)
                
    if (val$verbose) { print(val) }

    return (val)
}

rescale.GGP = function(data.ggp, x.mid, y.mid, x.rad, y.rad) {
    if (is.na(x.mid) & !is.na(x.rad)) stop("Please specify X midpoint if specifying X range ")
    if (is.na(y.mid) & !is.na(y.rad)) stop("Please specify Y midpoint if specifying Y range ") 
    if (!is.na(x.mid) & is.na(x.rad)) stop("Please specify X range if specifying X midpoint ")
    if (!is.na(y.mid) & is.na(y.rad)) stop("Please specify Y range if specifying Y midpoint ")

    if (is.na(x.mid) | is.na(y.mid)) {
        stop("Please define both x.mid and y.mid")
    }

    x.range.start = x.mid - x.rad
    x.range.end = x.mid + x.rad
    y.range.start = y.mid - y.rad
    y.range.end = y.mid + y.rad
    data.ggp = data.ggp + coord_cartesian(xlim=c(x.range.start, x.range.end), ylim=c(y.range.start, y.range.end))

    # in many cases it is much better to use coord_cartesian than xlim() like below - the latter cuts off data, e.g., geom_segment, which
    # is outside of range

#    if (!is.na(x.mid)) {
#        x.range.start = x.mid - x.rad
#        x.range.end = x.mid + x.rad
#        #data.ggp = data.ggp + xlim(range.start, range.end)  - this was the original call
#        data.ggp = data.ggp + coord_cartesian(xlim=c(x.range.start, x.range.end))
#    }
#    if (!is.na(y.mid)) {
#        y.range.start = y.mid - y.rad
#        y.range.end = y.mid + y.rad
#        #data.ggp = data.ggp + ylim(range.start, range.end)  - this was the original call
#        data.ggp = data.ggp + coord_cartesian(ylim=c(y.range.start, y.range.end))
#    }

    return(data.ggp)
}

args = parse_args()
if (!file.exists(args$ggp.fn)) {
    stop(paste("GGP file", args$ggp.fn, "not found."))
}
data.ggp = readRDS(args$ggp.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/

if (!is.na(args$special.x)) {
    data.ggp = data.ggp + geom_point(aes(x=args$special.x, y=args$special.y), color="red", alpha=0.75, shape=3, size=4, show_guide=FALSE)
}

if (args$no.flip) {
    data.ggp = rescale.GGP(data.ggp, args$x.mid, args$y.mid, args$x.rad, args$y.rad)
} else {
    data.ggp = rescale.GGP(data.ggp, args$y.mid, args$x.mid, args$y.rad, args$x.rad)
}

# set aspect ratio to be square.  Using theme() because coord_fixed() undoes coord_flip()
# not sure this does anything anyway
data.ggp = data.ggp + theme(aspect.ratio=1)

# Save
# By default R outputs PDF to Rplots.pdf.  To avoid littering we delete this.  Do not delete Rplots.pdf if that file
# already exists, though.
delete.Rplots = TRUE
if (file.exists("Rplots.pdf")) delete.Rplots=FALSE

# A dirty trick: ggsave will not save grobs, so mess with the the code directly.
# from http://stackoverflow.com/questions/18406991/saving-a-graph-with-ggsave-after-using-ggplot-build-and-ggplot-gtable
# ggsave = ggplot2::ggsave; body(ggsave) = body(ggplot2::ggsave)[-2]  

ggsave(plot=data.ggp, filename=args$out.fn, height=args$height, width=args$width, useDingbats=FALSE)
cat(sprintf("Saved to %s\n", args$out.fn))

if (delete.Rplots) {
    unlink("Rplots.pdf")
}

