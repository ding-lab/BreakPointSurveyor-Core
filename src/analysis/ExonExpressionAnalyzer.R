# Matthew Wyczalkowski and Jennifer Flynn
# m.wyczalkowski@wustl.edu
# version 3.0 - per gene analysis with permutation test across all genes of interest.  
#               Algorithm described in L_Expression/AlgorithmDetails.md of BPS.TCGA_Virus.Lite
#
# Usage: Rscript ExonExpressionAnalyzer.R [-v] [-S] [-s gene.out.fn] [-z zero.cutoff] [-n name] [-A nA] [-F] case.sample.id RPKM.dat exons.dat 
#
# Calculate and write to stdout p-value associated with gene expression in vicinity of integration event.
# Intra-integration event exons are retained.  Exons with excessive number of case/control samples where RSEM = 0 (50% by default)
# are discarded from further analysis.  Also calculate FDR, assuming all genes are in same test

# exons.dat: six-column BED file listing neighbor exons.  5th column must be 'downstream', 'upstream', or 'intra'
# RPKM.dat: RPKM data file with exons as rows and samples as columns
# sample.id: Sample ID of case, (format: TCGA-2F-A9KO-01)

# -v verbose output
# -s gene.out.fn: write gene p-values to TSV file
# -S: print genes data to stdout
# -z zero.cutoff: Discard exons for which "zero.cutoff" or more samples have RSEM value of zero [0.50]
# -n name: Optional integration event name used in output report [default is case.sample.id]
# -A nA: number of times to draw observations at random to create null distribution [10,000]
# -F: calculate FDR

library("reshape2")
library("plyr")

options("width"=170) # useful for debugging
options(warn=2)  # quit on warning


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

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    gene.out.fn = get_val_arg(args, "-s", NULL)
    print.genes = get_bool_arg(args, "-S")
    zero.cutoff = as.numeric(get_val_arg(args, "-z", 0.5))
    event.name = get_val_arg(args, "-n", NULL)
    nA = as.numeric(get_val_arg(args, "-A", 10000))
    do.FDR = get_bool_arg(args, "-F")

    # in reverse order
    exons.fn = args[length(args)];           args = args[-length(args)]
    RPKM.fn = args[length(args)];           args = args[-length(args)]
    sample.id = args[length(args)];           args = args[-length(args)]

    if (is.null(event.name)) {
        event.name = sample.id
    }

    val = list( 'verbose'=verbose, 'exons.fn'=exons.fn, 'RPKM.fn'=RPKM.fn, 'sample.id'=sample.id, 'gene.out.fn'=gene.out.fn,
                'zero.cutoff'=zero.cutoff, 'event.name'=event.name, 'print.genes'=print.genes, 'nA'=nA, 'do.FDR'=do.FDR)
    if (val$verbose) { print(val) }

    return (val)
}
# *** /Command Line Args ***

read.exons = function(exons.fn) {
    # Retain just the exons listed in exons.dat
    exons = read.table(exons.fn, col.names=c("chrom", "start", "end", "gene", "stream", "strand"))
    #  chrom    start      end gene   stream strand
    #  1    14 68086514 68086805 ARG2 upstream      +
    exons$exon.id = as.numeric(rownames(exons))
    return(exons)
}

# Read RPKM data, retain only entries of interest (as defined by exons), melt to long format, and indicate case/control
# Exons with an excessive number of zero-RPKM samples (as defined by zero.cutoff) are discarded.
# Format of returned RPKM data frame:

# Input RPKM data has the following four initial columns: chrom, start, end, name, with RPKM data in the following columns.
#   Header lines are mandatory.

# TODO
#   * Detail format out output data    

get.RPKM = function(RPKM.fn, exons.fn, zero.cutoff) {
    exons = read.exons(exons.fn)
    # Parse exon expression data
    RPKM.dat = read.table(RPKM.fn, row.names=NULL, header=TRUE, sep = "\t", check.names=FALSE, comment.char="#")
    if (nrow(RPKM.dat) == 0) {
        stop(paste("File", args$RPKM.fn, "has no data"))
    }
    RPKM.dat = RPKM.dat[RPKM.dat$chrom %in% exons$chrom & RPKM.dat$start %in% exons$start & RPKM.dat$end %in% exons$end,]

    # Add stream (upstream, downstream), strand, as well as exon ID information to RPKM.dat from exons
    RPKM.dat = merge(RPKM.dat, exons[,c("chrom", "start", "end", "gene", "stream", "strand", "exon.id")], by=c("chrom", "start", "end", "gene"))

    # chrom    start      end    gene strand TCGA-BA-4074-01A-01R-1436-07 TCGA-BA-4075-01A-01R-1436-07 TCGA-BA-4076-01A-01R-1436-07 
    # 14 68168602 68168652   RDH12      +                   0.00000000                   0.00000000                     0.000000
    RPKM = melt(RPKM.dat, id.vars=c("chrom", "start", "end", "gene", "stream", "strand", "exon.id"), variable.name = "barcode", value.name = "RPKM")
    RPKM$sample.id = substr(RPKM$barcode, 1, 15)
    RPKM$case = RPKM$sample.id == args$sample.id    # separate case from control
    if (!any(RPKM$case)) {
        stop(sprintf("Case %s has no RPKM data in %s.  Quitting.", args$sample.id, args$RPKM.fn))
    }
    if (length(unique(RPKM[RPKM$case,]$barcode)) > 1) {
        warning(sprintf("Case %s has more than one barcode.  Continuing.", args$sample.id))
        # If this happens will need to reevaluate how to calculate exon p-value
    }

    # Discard exons for which too many samples have RPKM=0.
    # Calculate fraction of samples with 0 expression for each exon and discard those which don't meet threshold
    RPKM.zeros = ddply(RPKM, "exon.id", summarize, zero.fraction = sum(RPKM==0)/length(RPKM) )
    exons = merge(exons, RPKM.zeros, by="exon.id", sort=FALSE)
    
    exons = exons[exons$zero.fraction < args$zero.cutoff,]
    if (nrow(exons) == 0) {
        stop(sprintf("No exons for sample %s after filtering.  Quitting.", args$sample.id))
    }

    RPKM = RPKM[RPKM$exon.id %in% exons$exon.id,]

    return(RPKM)
}

# difference.of.means is the the difference of means between x[A] and x[B].
# given vector x and index (or indices) A, calculate mean(x[A]) - mean(x[B]),
# where x[B] = x[-A] is x with indices A removed.
difference.of.means = function(A, x) {
    return(mean(x[A]) - mean(x[-A]))
}

# create random vector v of length l of elements [1..length(x)]
# evaluate abs(mean(x[v]) - mean(x[-v]))
random.abs.difference.of.means = function(l,x) {
    v = sample(length(x), l)
    return(abs(difference.of.means(v, x)))
}

# for every gene, calculate the P.value like,
#   p.value = count(delta[i] >= delta*)/nA
# where delta[i] = <X[Ai]> - <X[Bi]>
# with Ai, Bi the indices of randomly drawn case, control observations, i = 1..nA
# and delta* corresponding to real case, control 

get.gene.stats = function(exons, nA) {
# passed here are slices of RPKM.df which have the same gene name and stream, e.g.,
# chrom    start         end      gene     stream     strand exon.id    barcode                   RPKM       sample.id      case
#     9      14 68758600 68758697 RAD51B downstream      +       9 TCGA-BA-4074-01A-01R-1436-07  3.84690213 TCGA-BA-4074-01 FALSE

    # calculate span of gene from max, min positions of all exons.  Note that this will not span the integration event.
    gene.start = min(exons$start)
    gene.end = max(exons$end)

    number.exons=sum(exons$case)

    # We define delta as abs difference of average expression of exons in group A and B, where A is vector of randomly
    # chosen observations
    # we construct the null distribution of delta by drawing A observations at random nA times.
    null.deltas = replicate(nA, random.abs.difference.of.means(l=number.exons, x=exons$RPKM))

    # and define delta for the actual case
    case.delta = abs(difference.of.means(which(exons$case), exons$RPKM))
    case.epsilon = difference.of.means(which(exons$case), exons$RPKM)

    # We calculate p.value as the mean of p.value(>) and p.value(>=) to deal with case where
    # case RSEM = 0.  see notes 4/8/15
    c.gt = sum(null.deltas > case.delta)
    c.ge = sum(null.deltas >= case.delta)

    p.value = 0.5 * (c.gt + c.ge) / nA
    up.regulated = case.epsilon > 0

    return(data.frame( start = gene.start, end = gene.end, p.value = p.value, up.regulated = up.regulated))
}

# convenience function to reorder columns. From http://stackoverflow.com/questions/18339370/reordering-columns-in-a-large-dataframe
movetofirst = function(data, move) {
  data[c(move, setdiff(names(data), move))]
}

# write dataframe to disk, optionally append additional header line 
write.df = function(data, out.fn, header=NULL) {
    con = file(out.fn, open="wt")
    l = paste0("# Created ", Sys.time())
    if (!is.null(header)) {
        l = paste(header, l, sep="\n")
    }
    writeLines(l, con)
    write.table(data, con, sep="\t", quote=FALSE, row.names=FALSE)
    close(con)
    write( sprintf("    Saved to %s\n", out.fn), stderr())
}


args = parse_args()

RPKM.df = get.RPKM(args$RPKM.fn, args$exons.fn, args$zero.cutoff)
#  chrom    start      end gene   stream strand exon.id                      barcode      RPKM       sample.id  case
#  1    14 68086514 68086805 ARG2 upstream      +       1 TCGA-BA-4074-01A-01R-1436-07  9.041864 TCGA-BA-4074-01 FALSE

gene.Pvalue = ddply(RPKM.df, c("gene", "stream", "strand", "chrom"), function(exons) return(get.gene.stats(exons, args$nA)))
gene.Pvalue$integration.event = args$event.name
gene.Pvalue = movetofirst(gene.Pvalue, "integration.event")

# FDR calculations based on ~/src/genome/lib/perl/Genome/Model/Tools/Music/ClinicalCorrelation.pm.R
# method: Benjamini, Y., and Hochberg, Y. (1995). J Royal Stat Soc B 57, 289–300
if (args$do.FDR) {
    gene.Pvalue$FDR = p.adjust(gene.Pvalue$p.value, method="fdr") 
}

if (args$print.genes) {
    print(gene.Pvalue)
}
if (!is.null(args$gene.out.fn)) {
    h = sprintf("# Gene expression analysis for %s (%s)", args$sample.id, args$exons.fn)
    write.df(gene.Pvalue, args$gene.out.fn, h)
}


