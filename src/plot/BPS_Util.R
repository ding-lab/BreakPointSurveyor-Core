# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Common BreakpointSurveyor plot-related utilities.

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

# Parse chromosome region string.  Accepted formats: "-c 14", "-c 14:12345-12456"
# returns list ['range.pos', 'range.chr'], where range.pos is a list (e.g., [12345, 12456]) 
# and range.chr is a string ("14").  range.pos is NULL if not specified.
parse.range.str = function(chrarg) {
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


# save ggp object to either GGP binary representation or PDF file.
# if writing to PDF (pdf.out=TRUE), append ".pdf" to filename out.fn if it does not already have that extension
write.GGP = function(ggp, out.fn, pdf.out=FALSE) {
    # http://stat.ethz.ch/R-manual/R-devel/library/base/html/save.html
    if (!pdf.out) {
        cat(paste("Saving to GPP file", out.fn, "\n"))
        saveRDS(ggp, file=out.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
    } else {

        ext = file_ext(out.fn)
        if (ext != "pdf")
           out.fn = paste(sep=".", out.fn, "pdf") 

        cat(paste("Saving to PDF file", out.fn, "\n"))

        ggsave(plot=ggp, filename=out.fn, useDingbats=FALSE)
        unlink("Rplots.pdf")
    }
}
