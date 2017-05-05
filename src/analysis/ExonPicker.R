# Matthew Wyczalkowski and Jennifer Flynn
# m.wyczalkowski@wustl.edu
#
# Usage: Rscript ExonPicker.R [-v] [-K num.genes] [-E] -a ie.start -b ie.end -c ie.chrom -e exons.bed -o out.dat 

# Evaluate exons from genes upstream and downstream of integration event and write BED file describing these.
# We may select K such genes or all of them.
# 
# Mandatory parameters:
# -a, -b, -c: start, end, chrom of integration event.
# -e: exons BED file
# -o: output BED file (columns: chrom, start, end, gene.name, stream, strand)
#     stream is "upstream", "downstream", "intra" (inside of integration event)
# Optional parameters:
# -K: number of upstream, downstream genes to retain [default is to print all genes]
# -E: whether to include "extra" gene if integration event falls inside of a gene [default false]

#library("reshape2")
library("plyr")

options("width"=175) # useful for debugging

# *** Command Line Args ***
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
# Usage: Rscript ExonPicker.R [-v] [-K num.genes] [-E] -a ie.start -b ie.end -c ie.chrom -e exons.bed -o out.dat 
    verbose = get_bool_arg(args, "-v")
    E.bool = get_bool_arg(args, "-E")
    K = as.numeric(get_val_arg(args, "-K", NA))
    ie.start = as.numeric(get_val_arg(args, "-a", NA))
    ie.end = as.numeric(get_val_arg(args, "-b", NA))
    ie.chrom = get_val_arg(args, "-c", NULL)
    exons.bed.fn = get_val_arg(args, "-e", NULL)
    out.fn = get_val_arg(args, "-o", NULL)

    if (is.na(ie.start) | is.na(ie.end) | is.null(ie.chrom) )
        stop("Please pass all parameters -a ie.start -b ie.end -c ie.chrom")
    if (is.null(exons.bed.fn))
        stop("Please provide exons.bed filename")
    if (is.null(out.fn))
        stop("Please provide output filename")

    E = ifelse(E.bool, 1, 0)

    val = list( 'verbose' = verbose, 'E' = E, 'K' = K, 'ie.start' = ie.start, 
                'ie.end' = ie.end, 'ie.chrom' = ie.chrom, 'exons.bed.fn' = exons.bed.fn, 'out.fn' = out.fn)
    if (val$verbose) { print(val) }
    return (val)
}
# *** /Command Line Args ***

# be able to read bed files of arbitrary number of columns.  Columns named as described here:
#   http://bedtools.readthedocs.org/en/latest/content/general-usage.html
read.bed = function(bed.fn) {
    #data = read.csv(bed.fn, header = FALSE, sep = "\t", col.names=c("chrom", "start", "end", "name"), row.names=NULL)
    data = read.csv(bed.fn, header = FALSE, sep = "\t", row.names=NULL)
    names(data)[1:3] = c("chrom", "start", "end")
    if (length(names(data)) > 3) names(data)[4] = "name"
    if (length(names(data)) > 4) names(data)[5] = "score"
    if (length(names(data)) > 5) names(data)[6] = "strand"
    return(data)
}

# Read exons.  Optionally remove chrom not of interest
get.exon.bed = function(exons.bed.fn, chrom = NULL) {
    exons = read.bed(exons.bed.fn)
    if (!is.null(chrom))
        exons = exons[exons$chrom == chrom,]
    if (nrow(exons) == 0)
        return(NULL)
    #  chrom start   end    name score strand
    #  1     1 11868 12227 DDX11L1     .      +
    # Define exon pos, which is midpoint of each exon
    exons$pos = exons$start + (exons$start - exons$end)/2
    # rename column name, initialize $stream
    names(exons)[4] = "gene"
    exons$score = NULL
    exons$distance = NA
    return(exons)
}


get.exon.distance = function(exons, ie.start, ie.end) {
    # Define exon distance from integration event.  All exon positions evaluated at $pos
    # Distance calculated from exon midpoint to closest integration event boundary
    # exons lying within integration event have distance 0
    # downstream genes have distance > 0, upstream genes have distance < 0
    # downstream is defined in direction of increasing genomic position for + strand genes,
    #   and decreasing genomic position for - strand.  Upstream is the opposite.

    # Intra.event exons have distance 0
    # cases defined in doc/ExonExpressionSelectionWorkflow.pdf
#    e1 = exons[exons$pos >= ie.start & exons$pos <= ie.end,]
    intra.event = which(exons$pos >= ie.start & exons$pos <= ie.end)
    exons.A = which(exons$strand == "+" & exons$pos > ie.end)   # +'ve strand downstream
    exons.B = which(exons$strand == "+" & exons$pos < ie.start) # +'ve strand upstream
    exons.E = which(exons$strand == "-" & exons$pos > ie.end)   # -'ve strand upstream
    exons.D = which(exons$strand == "-" & exons$pos < ie.start) # -'ve strand downstream

    if (length(intra.event) != 0) exons[intra.event,"distance"] = 0
    if (length(exons.A) != 0) exons[exons.A,]$distance = exons[exons.A,]$pos - ie.end
    if (length(exons.B) != 0) exons[exons.B,]$distance = exons[exons.B,]$pos - ie.start
    if (length(exons.E) != 0) exons[exons.E,]$distance = ie.end - exons[exons.E,]$pos 
    if (length(exons.D) != 0) exons[exons.D,]$distance = ie.start - exons[exons.D,]$pos 
    return(exons)
}

# Select exons from K upstream and downstream genes nearest to exon 
# Exons have "stream" value of upstream, downstream, intra, distant
# Indicate which (if any) genes have exons or introns affected by integration event.
get.neighbors = function(exons, K, E, print.genes.affected=TRUE, skip.intra=FALSE) {
    exons$stream = "distant"
    intra.event = which(exons$dist == 0)
    if (length(intra.event) != 0) {
        exons[intra.event,]$stream = "intra"
        genes.exon.affected = unique(exons[intra.event,]$gene)
    } else {
        genes.exon.affected = NULL
    }

    # Find position of downstream genes, with gene position given by closest exon to integration event
    # we exclude intra-integration event exons so that genes which are entirely within integration event
    # do not count toward K genes outside of it.
    downstream.genes = ddply(exons[exons$distance > 0,], "gene", summarise, exon.dist = min(distance))
    upstream.genes = ddply(exons[exons$distance < 0,], "gene", summarise, exon.dist = min(-distance))

    # genes on both upsteam and downstream lists imply intragenic integration event.  If so, increase K by E
    # also evaluate whether genes have introns - but no exons - within integration event
    if (any(downstream.genes$gene %in% upstream.genes$gene)) {  # could also do intersect()
        if (!is.na(K))
            K = K + E
        genes.intron.affected = intersect(downstream.genes$gene, upstream.genes$gene)
        genes.intron.affected = setdiff(genes.intron.affected, genes.exon.affected)
    } else {
        genes.intron.affected = NULL
    }
    if (print.genes.affected) {
        cat(paste("Genes with exon in integration event:", paste(genes.exon.affected, collapse=" "), "\n"))
        cat(paste("Genes with intronic integration event:", paste(genes.intron.affected, collapse=" "), "\n"))
    }

    # Sort gene lists by distance and select the K closest ones
    # extract upstream and downstream exons, but not intra ones
    if (nrow(downstream.genes) > 0) {
        #downstream.genes = downstream.genes[order(downstream.genes$exon.dist),][1:K,]
        downstream.genes = downstream.genes[order(downstream.genes$exon.dist),]
        if (!is.na(K))
            downstream.genes = downstream.genes[1:K,]
        exons[exons$gene %in% downstream.genes$gene & exons$distance > 0,]$stream = "downstream"
    } 
    if (nrow(upstream.genes) > 0) {
        #upstream.genes = upstream.genes[order(upstream.genes$exon.dist),][1:K,]
        upstream.genes = upstream.genes[order(upstream.genes$exon.dist),]
        if (!is.na(K))
            upstream.genes = upstream.genes[1:K,]
        exons[exons$gene %in% upstream.genes$gene & exons$distance < 0,]$stream = "upstream"
    }
    return(exons)
}


args = parse_args()

cat(sprintf("integration event %d - %d\n", args$ie.start, args$ie.end))

exons = get.exon.bed(args$exons.bed.fn, args$ie.chrom) 
if (is.null(exons)) 
    stop(sprintf("No exons listed for chrom %s.  Exiting", args$ie.chrom))
    
exons = get.exon.distance(exons, args$ie.start, args$ie.end)
exons = get.neighbors(exons, args$K, args$E)
#      chrom    start      end          gene strand      pos  distance     stream
#      90451    14 64118014 64118217            U3      - 64117912 4515151.5 downstream
exons = exons[exons$stream != "distant",]

write.table(exons[,c("chrom", "start", "end", "gene", "stream", "strand")], file=args$out.fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat(sprintf("    Saved to %s\n", args$out.fn), file=stderr())



