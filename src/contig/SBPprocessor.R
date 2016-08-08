# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript SBPprocessor.R [-v] [-q qSBP.fn] [-o rSBP.fn] SBP.fn 
#
# Version 1.5 7/25/16 - Losing explicit reference to human and virus, using chromA and chromB instead
# Version 1.1 5/15/15 - outputs qSBP to indicate paired breakpoint information
# Version 1.0 5/12/15  
#
# process SBP (SAM BreakPoint) data to extract breakpoint positions as chrom A and B coordinates,
# where A < B by string comparison 
# Two output datasets are produced:
#   rSBP lists all unique breakpoints as chrom A, B coordinates
#   qSBP focuses on paired breakpoints (two or more breakpoints associated with one contig)
# 
# SBP.fn: List of SBP breakpoints.  
# -o rSBP.fn: output for reduced (R format) SBP data.  Writes to stdout by default
# -q qSBP.fn: output for paired breakpoints information.  Not generated if not specified.

# We define as "A first" those segment pairs where the first rname < second rname 
# This is appropriate for human/virus and human/human breakpoints and is consistent
# with ordering in BPC and BPR files.  This will not work reliably with intra-chromosomal
# translocations (where chrom A == chrom B)


options("width"=180) # useful for debugging
library(plyr)

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
    rSBP.fn = get_val_arg(args, "-o", "stdout")
    qSBP.fn = get_val_arg(args, "-q", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    SBP.fn = args[length(args)];      args = args[-length(args)]

    val = list('verbose'=verbose, 'rSBP.fn'=rSBP.fn,  'SBP.fn'=SBP.fn, 'qSBP.fn'=qSBP.fn)
    if (val$verbose) { print(val) }

    return (val)
}

# what a sample row of sbp.fn looks like:
#  /Users/mwyczalk/bin/examine_row ../A_Data/origdata/SBP/TCGA-CV-5971-01A-11D-1681-02.SBP.dat 3
#     1  query_name  8^128865748^gi|9626069|ref|NC_001357.1|^2902^CTX^143^-+.Contig-15.-1.-9.7.3.2.11.10.4.8.6.5.14.16
#     2  contig.id   1
#     3  bp.id   1
#     4  ref_name.Sa gi|9626069|ref|NC_001357.1|
#     5  is_primary.Sa   FALSE
#     6  is_forward.Sa   FALSE
# x   7  left_gpos.Sa    1
#     8  right_gpos.Sa   430
# x   9  left_rpos.Sa    3308
#    10  right_rpos.Sa   2879
#    11  ref_name.Sb 8
#    12  is_primary.Sb   TRUE
#    13  is_forward.Sb   TRUE
#    14  left_gpos.Sb    441
# x  15  right_gpos.Sb   1316
#    16  left_rpos.Sb    128865729
# x  17  right_rpos.Sb   128866602

# Read SBP file and mark rows in which chrom A first (which is unrelated to Sa and Sb).  
# Convert all positions (gpos and rpos) to 1-index format
# Retain only columns which are informative about breakpoints:
# discard positions not corresponding to breakpoint (the junction between Sa and Sb) and is_primary.* columns.
# query_name is also discarded here; it is referred to uniquely by contig.id.  
# finally, convert Sa/Sb quantities to A/B as appropriate
get.SBP.AB = function(sbp.fn) {

# A.first defines mapping between Sa/Sb and chrom A, B
# sbp.fn is in python-esque 0-index format.  Convert to R-ish 1-index format by adding 1 to all positions
    SBP = read.table(sbp.fn, sep="\t", header=TRUE)

    # construct sbp.AB, which has chrom A, B breakpoint positions
    # SBP$A.first = !grepl("^gi", SBP$ref_name.Sa)  - old, this is how virus was defined
    # note that typically chrom A will be human and chrom B virus.
    SBP$A.first = as.character(SBP$ref_name.Sa) < as.character(SBP$ref_name.Sb)

    drop.cols = c("left_gpos.Sa", "left_rpos.Sa", "right_gpos.Sb", "right_rpos.Sb", "query_name", "is_primary.Sa", "is_primary.Sb")
    SBP = SBP[,!(names(SBP) %in% drop.cols)]
    
    # convert remaining positions to 1-index format
    SBP$right_gpos.Sa = SBP$right_gpos.Sa + 1
    SBP$right_rpos.Sa = SBP$right_rpos.Sa + 1
    SBP$left_gpos.Sb = SBP$left_gpos.Sb + 1
    SBP$left_rpos.Sb = SBP$left_rpos.Sb + 1

    # now define A, B quantities based on whether chrom A or chrom B is Sa (requiring A < B by string comparison)
    SBP$A.rname = as.factor(ifelse(SBP$A.first, as.character(SBP$ref_name.Sa), as.character(SBP$ref_name.Sb)))
    SBP$B.rname = as.factor(ifelse(SBP$A.first, as.character(SBP$ref_name.Sb), as.character(SBP$ref_name.Sa)))
    SBP$A.bp.rpos = ifelse(SBP$A.first, SBP$right_rpos.Sa, SBP$left_rpos.Sb)
    SBP$B.bp.rpos = ifelse(SBP$A.first, SBP$left_rpos.Sb, SBP$right_rpos.Sa)
    SBP$A.bp.gpos = ifelse(SBP$A.first, SBP$right_gpos.Sa, SBP$left_gpos.Sb)
    SBP$B.bp.gpos = ifelse(SBP$A.first, SBP$left_gpos.Sb, SBP$right_gpos.Sa)
    SBP$A.is_forward = ifelse(SBP$A.first, SBP$is_forward.Sa, SBP$is_forward.Sb)
    SBP$B.is_forward = ifelse(SBP$A.first, SBP$is_forward.Sb, SBP$is_forward.Sa)

    # and drop unneeded quantities
    drop.cols = c("ref_name.Sa", "is_forward.Sa", "right_gpos.Sa", "right_rpos.Sa", "ref_name.Sb", "is_forward.Sb", "left_gpos.Sb", "left_rpos.Sb")
    SBP = SBP[,!(names(SBP) %in% drop.cols)]

    return(SBP)
}

# process SBP.AB data to obtain the reduced rSBP format.  
# and only unique breakpoints - as defined by their reference position pairs - are retained
get.rSBP = function(SBP.AB) {
    return(unique(SBP.AB[,c("A.rname", "A.bp.rpos", "B.rname", "B.bp.rpos")]))
}

# See BreakPointParser.R, mark.segments() for details of how this works - the implementation of combining
# paired breakpoints is similar to combining segment pairs 
mark.pairs = function(SBP.AB) {
    num.SBP = nrow(SBP.AB)
    SBP.AB = SBP.AB[with(SBP.AB, order(bp.id)), ]

    # define columns which will help combine individual breakpoints into pairs of left (bpM), right (bpN) breakpoints 
    SBP.AB$pair.id = seq(0,num.SBP-1)

    # stuff here deals with expanding mid BP (those with another BP on each end) into two pairs
    SBP.AB$A.pos = rep.int("right", num.SBP)
    SBP.AB$seg.type = rep.int("mid", num.SBP)

    # now take care of first, last entries
    SBP.AB$pair.id[1] = 1
    SBP.AB$A.pos[1] = "left"
    SBP.AB$seg.type[1] = "first"
    SBP.AB$seg.type[num.SBP] = "last"

    # "mid" entries are duplicated to be on left in next bp
    SBP.mid = SBP.AB[SBP.AB$seg.type == "mid",]
    if (nrow(SBP.mid) > 0) {
        SBP.mid$A.pos = "left"
        SBP.mid$pair.id = SBP.mid$pair.id + 1
        SBP.AB = rbind(SBP.AB, SBP.mid)
    }

    return(SBP.AB)
}

# qSBP reports on breakpoint pairs, wherein two (or more) breakpoints occur on same contig.
#   for a contig with N > 1 breakpoints, N-1 qSBP rows generated
# In each breakpoint pair, breakpoints are labeled bpM and bpN.
get.qSBP = function(SBP) {

# returns:
# contig.id
# bpM.id, bpN.id - correspond to bp.id in SBP 
# bpM.A.rpos, bpM.B.rpos, bpN.A.rpos, bpN.B.rpos
# bpM.left.is_A, bpN.left.is_A - indicates whether left segment of each breakpoint is chrom A

    # get count of repeated contig.ids.  http://stackoverflow.com/questions/18201074/find-how-many-times-duplicated-rows-repeat-in-r-data-frame
    contig.count = ddply(SBP,"contig.id",nrow)
    # contig.id V1
    # 1  1
    # 6  2
    repeated.contig.ids = contig.count[contig.count$V1 > 1,"contig.id"]
    SBP.paired = SBP[SBP$contig.id %in% repeated.contig.ids,]
    if (nrow(SBP.paired) == 0) {
        return(NULL)
    }
#    contig.id bp.id A.first A.rname                       B.rname A.bp.rpos B.bp.rpos A.bp.gpos B.bp.gpos A.is_forward B.is_forward
#6           6     1        TRUE      14 gi|310698439|ref|NC_001526.2|  68741695      3434       104       103         TRUE         TRUE
#7           6     2       FALSE      14 gi|310698439|ref|NC_001526.2|  68683566      3979       645       648         TRUE         TRUE

    # Our aim now is to pair SBP (breakpoint) lines to create qSBP breakpoint pairs.  This approach is very similar to that in
    # parse.pSBP in BreakPointParser.R: each breakpoint in one contig is assigned a pair.id, a horizontal position (bpM or bpN), and 
    # a bp.type (first, mid, or last).  "mid" breakpoints occur in contigs with three or more breakpoints, and need to be paired 
    # with a preceeding breakpoint as well as a following breakpoint to form two pairs; to do this, they are duplicated.

    SBP.paired = ddply(SBP.paired, "contig.id", function(df) return( mark.pairs(df) ) )

    # now split into left (Sa) and right (Sb) segments, clean up, and merge
    SBP.bpM = SBP.paired[SBP.paired$A.pos=="left",]  
    SBP.bpN = SBP.paired[SBP.paired$A.pos=="right",]   # right segment

    drop.cols = c("A.bp.gpos", "B.bp.gpos", "seg.type", "A.pos")
    SBP.bpM = SBP.bpM[,!(names(SBP.bpM) %in% drop.cols)]
    SBP.bpN = SBP.bpN[,!(names(SBP.bpN) %in% drop.cols)]

    qSBP = merge(SBP.bpM, SBP.bpN, by=c("contig.id", "pair.id"), suffixes=c(".bpM", ".bpN"))

    return(qSBP)
}


args = parse_args()

SBP.AB = get.SBP.AB(args$SBP.fn)
rSBP = get.rSBP(SBP.AB)
qSBP = get.qSBP(SBP.AB)

# write rSBP to either file or stdout
if (args$rSBP.fn == "stdout") {
    f = stdout()
} else {
    cat(paste("Writing rSBP to", args$rSBP.fn, "\n"))
    f = file(args$rSBP.fn)
}

# with col.names removed, this is in BPC format
write.table(rSBP, file = f, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# write qSBP only if interested in it.
if (!is.null(args$qSBP.fn)) {
    if (is.null(qSBP)) {
        cat(paste("No paired breakpoints for ", args$qSBP, "\n"))
    } else {
        cat(paste("        qSBP to", args$qSBP.fn, "\n"))
        write.table(qSBP, file = args$qSBP.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
}

