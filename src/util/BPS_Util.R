# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Common BreakpointSurveyor utilities.

library("tools")

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

# split a string like "123-456" into a numeric list [123,456]
get.range = function(range.str) {
    range.pos=as.numeric(strsplit(range.str, "-")[[1]])
    if (any(is.na(range.pos))) {
        stop("Error parsing: ", range.str)
    }
    return(range.pos)
}

# Parse chromosome region string.  Accepted formats: NULL, "14", "14:12345-12456", "12345-12456"
# returns list ['range.pos', 'range.chr'], where range.pos is a list (e.g., [12345, 12456]) 
# and range.chr is a string ("14").  range.chr and range.pos are NULL if not specified.
parse.range.str = function(chrarg) {
    range.pos = NULL
    range.chr = NULL
    if (!is.null(chrarg)) {
        if (grepl(":", chrarg)) {  # C:A-B
            chrlist = strsplit(chrarg, ":")[[1]]
            range.chr = chrlist[1]
            range.pos = get.range(chrlist[2])
        } else if (grepl("-", chrarg)) {  # A-B
            range.pos = get.range(chrarg)  # C
        } else {
            range.chr = chrarg
        }
    } 
    return( list('range.pos'=range.pos, 'range.chr'=range.chr) )
}

# reads in a BPC file which specifies breakpoint coordinates with optional attributes column.
# flip.ab is an optional attribute, if true swaps columns A <-> B
read.BPC = function(BPC.fn, flip.ab=FALSE) {
    data = read.table(BPC.fn, sep='\t')
    if (ncol(data) == 4) {
        if (!flip.ab) {
            names(data) = c("chrom.A", "pos.A", "chrom.B", "pos.B")
        } else {
            names(data) = c("chrom.B", "pos.B", "chrom.A", "pos.A")
        }
        data$attribute = "default"
    } else if (ncol(data) == 5) {
        if (!flip.ab) {
            names(data) = c("chrom.A", "pos.A", "chrom.B", "pos.B", "attribute")
        } else {
            names(data) = c("chrom.B", "pos.B", "chrom.A", "pos.A", "attribute")
        }
    } else {
        stop("Unexpected number of rows in ", BPC.fn)
    }
    data$attribute = factor(data$attribute)
    return(data)
}

# * BPR: chrom.A, pos.A.start, pos.A.end, chrom.B, pos.B.start, pos.B.end, [attribute] 
# flip.ab is an optional attribute, if true swaps columns A <-> B
read.BPR = function(BPR.fn, flip.ab = FALSE) {
    data = read.table(BPR.fn, sep='\t')
    if (ncol(data) == 6) {
        if (!flip.ab) {
            names(data) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end")
        } else {
            names(data) = c("chrom.B", "pos.B.start", "pos.B.end", "chrom.A", "pos.A.start", "pos.A.end")
        }
        data$attribute = "default"
    } else if (ncol(data) == 7) {
        if (!flip.ab) {
            names(data) = c("chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "attribute")
        } else {
            names(data) = c("chrom.B", "pos.B.start", "pos.B.end", "chrom.A", "pos.A.start", "pos.A.end", "attribute")
        }
    } else {
        stop("Unexpected number of rows in ", BPR.fn)
    }
    data$attribute = factor(data$attribute)
    return(data)
}

filter.BPC = function(data, range.chr.A, range.pos.A, range.chr.B=NULL, range.pos.B=NULL) {
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

# be able to read bed files of arbitrary number of columns.  Columns named as described here:
#   http://bedtools.readthedocs.org/en/latest/content/general-usage.html
read.BED = function(bed.fn) {
    # this complains if file has zero lines
    data = read.csv(bed.fn, header = FALSE, sep = "\t", row.names=NULL)
    names(data)[1:3] = c("chrom", "start", "end")
    if (length(names(data)) > 3) names(data)[4] = "name" 
    if (length(names(data)) > 4) names(data)[5] = "score" 
    if (length(names(data)) > 5) names(data)[6] = "strand" 
    return(data)
}

filter.BED = function(data, range.chr, range.pos) {
    # keep rows which overlap with range of interest.
    if (!is.null(range.chr)) data = data[data$chrom %in% range.chr,]
    if (!is.null(range.pos)) {
        in.range = is.overlapping(data$start, data$end, range.pos[1], range.pos[2])
        data = data[in.range,]
    } 
    return(data)
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

