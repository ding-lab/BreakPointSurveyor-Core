# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript RPKM_Joiner.R [-v] barcodes.dat out.fn

# Process RPKM files given by barcodes.dat and combine column-wise into one data file
# barcodes.dat is a filename with list of all barcodes and their data files
# Exon information is retained, RNA expression values are concatenated column-wise, and data file is
# written to out.fn
# chromosome names have 'chr' removed, i.e., 'chr2' becomes '2'.  This is required for consistency
# with BED files in this workflow.
# older versions (with 'chr') can be fixed with, `sed 's/^chr//' old.dat > new.dat`
#   (though this does clobber the header)

library("reshape2")
library("plyr")
library("data.table")
#library("ggplot2")

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

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.fn = args[length(args)];               args = args[-length(args)]
    barcodes.fn = args[length(args)]; args = args[-length(args)]

    val = list( 'verbose'=verbose, 'barcodes.fn'=barcodes.fn, 'out.fn'=out.fn)
    if (val$verbose) { print(val) }

    return (val)
}

# parse RPKM data file and return data frame
# retain barcode
get_rpkm_data = function(fn, barcode) {
    data.bar = read.table(fn, col.names=c("chrom", "start", "end", "exon.code", "rpkm"))
    #  chrom     start       end          exon.code       rpkm
    #  1  chr2 163200847 163201035     GCA:upstream:+ 0.44133959
    data.bar$chrom = gsub("chr", "", data.bar$chrom, fixed=TRUE)  # remove "chr" from chromosome name
    data.bar$barcode = barcode
    # using sample.id is problematic because get multiple same values sometimes.  Barcode seems unique.
    #data.bar$sample.id = substr(barcode, 1, 15)  # e.g., TCGA-B7-5816-01
    data.bar.table=data.table(data.bar)   # this is the key to making this faster.  read.table could be fread...
    return(data.bar.table)
}

args = parse_args()

barcode.list = read.table(args$barcodes.fn, row.names=NULL, header=FALSE, col.names=c("name", "path"), sep="\t")
data = NULL

# Note: for large datasets (more several hundred or more) this becomes quite slow, because rbind
# is memory inefficient.  The way to do this is with data.table:
# http://stackoverflow.com/questions/11486369/growing-a-data-frame-in-a-memory-efficient-manner
# More about data.table: http://www.marketingdistillery.com/2014/05/05/working-with-large-data-sets-in-r-the-data-table-module/
# concatenation example here: http://zevross.com/blog/2013/11/25/data_table/
for (i in 1:nrow(barcode.list)) {
    fn = as.character(barcode.list[i,]$path)
    b = barcode.list[i,]$name
    cat(paste(b, "\n"))
    d = get_rpkm_data(fn, b)
    data = rbindlist( list(data, d)) 
}

# song's output fields are dependent on BED file input format.  Here is the current format:
#   chrom     start       end          exon.code       rpkm       sample.id
#   1:     2 163200847 163201035     GCA:upstream:+ 0.44133959 TCGA-2F-A9KO-01

# New format does not have additional exon code fields, only has gene.
#   chrom    start      end exon.code       rpkm                      barcode
#   1:    14 68026324 68026408   PLEKHH1  0.1708560 TCGA-BA-4074-01A-01R-1436-07


data.wide=dcast(unique(data), chrom + start + end + exon.code ~ barcode, value.var="rpkm")

# http://stackoverflow.com/questions/4350440/using-strsplit-with-data-frames-to-split-label-columns-into-multiple
exon.dets = do.call(rbind, strsplit(as.character(data.wide$exon.code), ":"))
data.wide$gene = exon.dets[,1]
#data.wide$strand = exon.dets[,3]  # ignoring strand of gene
data.wide$exon.code = NULL

# Rearrange columns for convenience
# from http://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe
movetofirst <- function(data, move) {
#  data[c(setdiff(move, names(data)), move)]
data[c(move, setdiff(names(data), move))]

}

# no longer have strand info
#data.wide = movetofirst(data.wide, c("chrom", "start", "end", "gene", "strand"))
data.wide = movetofirst(data.wide, c("chrom", "start", "end", "gene"))

con = file(args$out.fn, open="wt")
timestamp = paste0("# Created ", Sys.time())
writeLines(paste("# RPKM data ", timestamp, sep="\n"), con)
write.table(data.wide, con, sep="\t", quote=FALSE, row.names=FALSE)
close(con)
cat(sprintf("    Saved to %s\n", args$out.fn))
