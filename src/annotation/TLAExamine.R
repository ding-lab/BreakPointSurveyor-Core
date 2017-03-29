# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript TLAExamine.R [-v] [-m mode] [-f dat.fn] [-o out.fn] 
# Version 2.0 (4-27-15)

# Processes GTF and VCF files by expanding a column of key/value pairs into multiple columns with value field
# Values with repeated keys are separated by ';'
# -f dat.fn: an ENSEMBL GTF or TIGRA VCF file.  If not specified reads from stdin
# -o out.fn: TSV format output.  If not specified writes to stdout
# -m mode: type of file, either "vcf" or "gtf".  If not specified then based on dat.fn extension

suppressMessages(library(splitstackshape))
suppressMessages(library("reshape2"))
library(tools)

options("width"=180) # useful for debugging

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

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    out.fn = get_val_arg(args, "-o", NULL)
    dat.fn = get_val_arg(args, "-f", NULL)
    file.type = get_val_arg(args, "-m", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    # dat.fn = args[length(args)];               args = args[-length(args)]

    val = list( 'verbose'=verbose, 'out.fn'=out.fn, 'dat.fn'=dat.fn, 'file.type'=file.type)
    if (val$verbose) { print(val) }

    return (val)
}

trim = function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

# obtain and validate file format 
# if file.type is passed use that, otherwise obtain file format from extension of data.fn
# if format is unknown, or if data.fn extension is ambiguous, quit with an error
detect.file.format = function(file.type, data.fn) {
    if (is.null(file.type)) {
        file.type = file_ext(data.fn)
    }
    if (!file.type %in% c("vcf", "gtf")) {
        stop(paste("Unknown file type", file.type))
    }
    return(file.type)
}

cast.wide = function(data.split, file.format) {
    if (file.format == 'gtf') {
        data = dcast(data.split, seqname + source + feature + start + end + score + strand + frame + feature.id ~ key, value.var="value")
    } else if (file.format == 'vcf') {
        data = dcast(data.split, chrom + pos + id + ref + alt + qual + filter + format + spikein + V11 + feature.id ~ key, value.var="value")
    }
}

args = parse_args()

file.format = detect.file.format(args$file.type, args$dat.fn)

if (file.format == "gtf") {
    # col.names from ENSEMBL GTF specification, see D_Ensembl/README.ensembl
    all.cols = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    col.complex = "attribute"
    field.delimit = ";"
    key.value.delimit = " "
} else if (file.format == "vcf") {
    # VCF colnames from vafFilter.py
    all.cols = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "spikein", "V11")
    col.complex = "info"
    field.delimit = ";"
    key.value.delimit = "="
}

if (is.null(args$dat.fn)) {
    fin = file('stdin')
} else {
    fin = args$dat.fn
}

data = read.table(fin, col.names=all.cols, sep="\t", row.names=NULL)
data$feature.id = row.names(data)  # to ensure uniqueness when doing dcast

# cSplit is the magic here
# http://www.r-bloggers.com/splitstackshape-v1-4-0-for-r/
# http://stackoverflow.com/questions/24595421/how-to-strsplit-data-frame-column-and-replicate-rows-accordingly
data.split = as.data.frame(cSplit(data, col.complex, sep = field.delimit, direction="long"))


key.value = as.data.frame(do.call(rbind, strsplit(trim(data.split[,col.complex]), key.value.delimit)))
names(key.value) = c("key", "value")
data.split = cbind(data.split, key.value)
data.split[,col.complex] = NULL

# http://stackoverflow.com/questions/14262741/combining-duplicated-rows-in-r-and-adding-new-column-containing-ids-of-duplicate
# collapse duplicate attribute entries into one row, separated by semicolons.  An example where this is needed is in a GTF attribute
# string like, 
#   tag "cds_end_NF"; tag "mRNA_end_NF";
data.split = aggregate(data.split, data.split[,c("feature.id", "key")], FUN = function(X) paste(unique(X), collapse=";"))
# aggregate duplicates two columns, so we get rid of them.
data.split = data.split[,3:ncol(data.split)]

# Cast into wide format
data = cast.wide(data.split, file.format)

# remove the feature ID column because its for internal use only
data$feature.id = NULL
data = unique(data)

if (is.null(args$out.fn)) {
    out.fn = ""
} else {
    out.fn = args$out.fn
}

write.table(data, file=, sep="\t", quote=FALSE, row.names=FALSE)
if (is.null(args$out.fn)) {
    cat(sprintf("    Saved to %s\n", args$out.fn))
}
