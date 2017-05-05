# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Usage: Rscript Pindel_RP.Reader.R [-v] [-V] [-p virus.prefix] [-S] PindelRP.dat out.BPR
#
#   Create Breakpoint Region file (BPR) based on output of Pindel RP module.
#
#   Optionally retain only inter-chromosomal translocations or human-virus breakpoints
#   Attribute (strand) column of BPR is string composed of chrom A and B strand info, e.g., "A+ B-"
#
# Input arguments:
#
# * PindelRP.dat: Output of Pindel RP module.  If 'stdin', read from stdin
# * out.BPR: BPR file representing Pindel RP breakpoint prediction regions
#
# Optional arguments:
#
# -v: Verbose output
# -S: Exclude intra-chromosomal events 
# -V: Retain only human-virus breakpoints, with virus "chromosome" identified by virus.prefix, below
# -p virus.prefix: string identifying virus reference.  Default 'gi'
#
# BPR and BPC file format descriptions: doc/FileFormat.txt

options("width"=180) # useful for debugging

source_relative = function(source.fn) {
    # http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
    initial.options = commandArgs(trailingOnly = FALSE)
    file.arg.name = "--file="
    script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename = dirname(script.name)
    other.name = paste(sep="/", script.basename, source.fn)
    source(other.name)
}
source_relative("../util/BPS_Util.R")

options("width"=300) # change terminal window width
#library("bitops")
#library(data.table)

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    exclude.intra = get_bool_arg(args, "-S")
    filter.virus = get_bool_arg(args, "-V")
    virus.prefix = get_val_arg(args, "-p", "gi")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.fn = args[length(args)]; args = args[-length(args)]
    dat.fn = args[length(args)]; args = args[-length(args)]

    val = list('verbose'= verbose, 'out.fn'= out.fn, 'dat.fn'=dat.fn, 'exclude.intra'=exclude.intra, 
            'filter.virus'=filter.virus, 'virus.prefix'=virus.prefix)
    if (val$verbose) { print(val) }
    return (val)
}

read.Pindel_RP = function(pindel.fn) {
    pindel.df = read.table(pindel.fn, col.names=c("chrom1", "start1", "end1", "strand1", "length1", 
            "chrom2", "start2", "end2", "strand2", "length2", 'event.size', "support"), sep="\t",
            colClasses = c("character", 'numeric', 'numeric', 'character', 'numeric', 'character', 
            'numeric', 'numeric', 'character', 'numeric', 'numeric', 'character'))
    if (nrow(pindel.df) == 0) 
        pindel.df = NULL
    return(pindel.df)
}

# optionally delete cases where chrom1 = chrom2
# optionally include only cases where one chrom starts with virus.prefix
filter.pindel = function(pindel.df, exclude.intra = FALSE, filter.virus = FALSE, virus.prefix="gi") {
    if (is.null(pindel.df)) return(NULL)

    if (exclude.intra)
        pindel.df = pindel.df[pindel.df$chrom1 != pindel.df$chrom2,]

    if (filter.virus) {
        pre = paste0("^", virus.prefix)
        keepA = which(grepl(pre,pindel.df$chrom1) & !grepl(pre,pindel.df$chrom2))
        keepB = which(!grepl(pre,pindel.df$chrom1) & grepl(pre,pindel.df$chrom2))
        pindel.df = pindel.df[ c(keepA, keepB), ]
    }

    if (nrow(pindel.df) == 0) {
        pindel.df = NULL
    }
    return(pindel.df)
}

as.BPR = function(pindel.df) {
    if (is.null(pindel.df)) return(NULL)
# * BPR: chrom.A, pos.A.start, pos.A.end, chrom.B, pos.B.start, pos.B.end, [attribute]
#  Require chrom.A < chrom.B (or, posA.start < posB.start if chrom.A = chrom.B)

    A.first = which(pindel.df$chrom1 < pindel.df$chrom2)
    A.first.pos = which( (pindel.df$chrom1 == pindel.df$chrom2) & (pindel.df$start1 < pindel.df$start2) )
    A.first = sort(unique(c(A.first, A.first.pos)))
    B.first = setdiff( seq_len(nrow(pindel.df)), A.first )

    pindel.BPR.A = pindel.df[A.first, c("chrom1", "start1", "end1", "strand1", "chrom2", "start2", "end2", "strand2")]
    pindel.BPR.B = pindel.df[B.first, c("chrom2", "start2", "end2", "strand2", "chrom1", "start1", "end1", "strand1")]
    names(pindel.BPR.A) = c("chrom.A", "pos.A.start", "pos.A.end", "strand.A", "chrom.B", "pos.B.start", "pos.B.end", "strand.B")
    names(pindel.BPR.B) = c("chrom.A", "pos.A.start", "pos.A.end", "strand.A", "chrom.B", "pos.B.start", "pos.B.end", "strand.B")
    pindel.BPR = rbind(pindel.BPR.A, pindel.BPR.B)

    # Create strand string
    pindel.BPR$strand = sprintf("A%s B%s", pindel.BPR$strand.A, pindel.BPR$strand.B)
    pindel.BPR$strand.A = NULL
    pindel.BPR$strand.B = NULL

    return(pindel.BPR)
}

args = parse_args()

# To read in from stdin.
f = if (args$dat.fn == "stdin") file("stdin") else args$dat.fn

#    f.in = file("stdin")
#    data = read.csv(f, sep='\t', header=FALSE, quote="")  # the " character can occur in SAM file and quoting must be disabled

pindel.df = read.Pindel_RP(f)
pindel.df = filter.pindel(pindel.df, args$exclude.intra, args$filter.virus, args$virus.prefix)
pindel.BPR = as.BPR(pindel.df)

### need to put a # character before header column
# For convenience and to identify the attribute column, write header as a comment line
cn = c("# chrom.A", "pos.A.start", "pos.A.end", "chrom.B", "pos.B.start", "pos.B.end", "strand")
write.table(pindel.BPR, file=args$out.fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=cn) #, row.names=FALSE)

# FYI: code to write arbitrary headers, including timestamp.  
#con = file(args$out.fn, open="wt")
#timestamp = paste0("# Created ", Sys.time())
#writeLines(paste("# Data File Description ", timestamp, sep="\n"), con)
#write.table(discovery.join, con, sep="\t", quote=FALSE, row.names=FALSE)
#close(con)
#cat(sprintf("    Saved to %s\n", args$out.fn))




