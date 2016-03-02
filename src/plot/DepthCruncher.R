# Usage: Rscript DepthCruncher.R [-v] [-C] [-P] [-s discordant.dat] [-p pindel_RP.dat] [-R HHRP] [-r rSBP.dat] [-V] [-L] [-u num.reads] [-l read.length]
#                                 -a chrom.range.start -b chrom.range.end -c chrom -n integration.event.name depth.dat depth.ggp
#
# visualize read depth of chrom or virus around one integration event.  Overlay Tigra breakpoints, discordant read pairs,
# and read depth(CBS) segmentation when available.  Save to GGP file.
# -C: skip CBS calculation and plotting (speeds up plotting for testing purposes)
# -P: save as PDF instead of GGP 
# -V: this is a virus (as opposed to chrom).  Relevant for reading discordant read, pindel, and tigra data
# -a chrom.range.start, -b chrom.range.end, -c chrom: defines chrom region of interest 
#    Typically this is the integration event plus flanking region, e.g. +/- 50Kbp.  
#    Values used to filter discordant reads, pindel and tigra data, as well as chrom depth.  Virus depth is not filtered.
# -n name: unique name of integration event (e.g., TCGA-BA-4077-blah.chr14A
# -s BPC.dat: optional processed BPC file with discordant reads 
#       breakpoint.BPC.fn is a breakpoint BPC file with chrom.A < chrom.B (alphabetical order).  
# -p pindel_RP.dat: pindel RP data file
# -R HHRP: HHRP file indicating human-human breakpoints
# -r rSBP: file indicating SBP-detected breakpoints
# -L: Plot read depth as log2(depth / <depth>), where <depth> is the median read depth
# -u: number of reads (typically mapped) in BAM file.  Necessary for normalizing to copy number
# -l: read length.  Necessary for normalizing to copy number

# Significant overlap with /Users/mwyczalk/Data/Virus/Virus_2013.9a/CombinedPlot/BreakpointCompare/B_BreakpointPlots.PBP/src/BreakpointPlotter.R
# version 2.0 6/4/15  Adds support for HHRP and SBP data.  Tigra support deleted.  
#                     Effort is made so homogeneize data structures and procedures with that in BreakpointPlotter v3.0.  In particular,
#                     do not coerce all data into one structure.  
# version 2.1 6/29/15 copy number is implemented

suppressPackageStartupMessages(library("ggplot2"))
options("width"=300) # change terminal window width
library("reshape2")

library(DNAcopy)  # DNAcopy is installed with source("http://bioconductor.org/biocLite.R"); biocLite(); biocLite("DNAcopy")

set.seed(25)
options(scipen=999)

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
    skip.cbs = get_bool_arg(args, "-C")
    is.virus = get_bool_arg(args, "-V")
    save.pdf = get_bool_arg(args, "-P")
    chrom.range.start = get_val_arg(args, "-a", NA)  # numeric quantities are NA so they can be tested with is.na()
    chrom.range.end = get_val_arg(args, "-b", NA)
    chrom = get_val_arg(args, "-c", NULL)
    integration.event.name = get_val_arg(args, "-n", NULL)
    breakpoint.BPC.fn = get_val_arg(args, "-s", NULL)
    pindel.fn = get_val_arg(args, "-p", NULL)
    HHRP.fn = get_val_arg(args, "-R", NULL)
    rSBP.fn = get_val_arg(args, '-r', NULL)
    plot.log.depth = get_bool_arg(args, "-L")
    num.reads = as.numeric(get_val_arg(args, "-u", NA))
    read.length = as.numeric(get_val_arg(args, "-l", NA))

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    ggp.fn = args[length(args)]; args = args[-length(args)]
    depth.fn = args[length(args)];      args = args[-length(args)]

    val = list( 'ggp.fn'=ggp.fn, 'HHRP.fn'=HHRP.fn, 'rSBP.fn'=rSBP.fn, 'depth.fn' = depth.fn, 'is.virus'=is.virus,
                'breakpoint.BPC.fn'=breakpoint.BPC.fn, 'verbose'=verbose, 'pindel.fn'=pindel.fn, 
                'skip.cbs'=skip.cbs, 'chrom.range.start'=as.numeric(chrom.range.start), 'chrom.range.end'=as.numeric(chrom.range.end),
                'chrom'= chrom, 'integration.event.name'= integration.event.name, 'save.pdf'=save.pdf,
                'plot.log.depth'=plot.log.depth, 'num.reads'=num.reads, 'read.length'=read.length)
    if (val$verbose) { print(val) }

    return (val)
}

calculateCBS = function(cn, do.log) {
# CBS code from /Users/mwyczalk/Data/Virus/Virus_2013.9a/ReadDepth/CBS/CBS.R
# Example of samtools depth output:
#14  68649196    6
#14  68649210    4
#14  68649224    8

    #CNA.object = CNA( log.depth = log2(cn$depth/median(cn$depth,na.rm=T)), chrom = cn$chrom, maploc = cn$pos, data.type = 'logratio')
    if (do.log) {
        CNA.object = CNA( genomdat = log2(cn$depth/median(cn$depth,na.rm=T)), chrom = cn$chrom, maploc = cn$pos)
    } else {
        print("CBS without log scaling is experimental")
        # Not certain this will work
        CNA.object = CNA( genomdat = cn$plot.depth, chrom = cn$chrom, maploc = cn$pos)
    }
    CNA.object = smooth.CNA(CNA.object)  
    d = segment(CNA.object, verbose=0, min.width=5, undo.splits="sdundo", undo.SD=3, alpha=1e-5)  # Smaller alpha makes fewer segments
    segs = d$output
    #write.table(segs[,2:6], file=outfn, row.names=F, col.names=F, quote=F, sep="\t")
    return (segs[,2:6])
}

# reads in a discordant read file as generated by A_Data/src/DiscordantCleaner.R
# Optionally filter by chrom, chrom position range, and virus name
# returns data frame with columns ("virus.name", "virus.pos", "chrom.name", "chrom.pos") 
#   Note this function differes significantly from function with same name in BreakpointPlotter.R v2.0

# This is old and is being replaced by read.BPC
read.discordant = function(breakpoint.BPC.fn, chrom=NULL, chrom.range.start=NA, chrom.range.end=NA, virname=NULL) {
    data = read.table(breakpoint.BPC.fn, sep='\t', header=TRUE, quote="") # the " character can occur in SAM file and quoting must be disabled

    # Specify chromosomes here.
    if (! is.null(chrom)) {
        data = data[data$chrom.name %in% chrom,]
    }
    if (!is.na(chrom.range.start) & !is.na(chrom.range.end)) {
        data = data[data$chrom.pos >= chrom.range.start & data$chrom.pos <= chrom.range.end,]
    } else if (!is.na(chrom.range.start)) {
        data = data[data$chrom.pos >= chrom.range.start,]
    } else if (!is.na(chrom.range.end)) {
        data = data[data$chrom.pos <= chrom.range.end,]
    }
    if (! is.null(virname) ) {
        #cat(paste("Keeping only virus", virname, "\n"))
        data = data[data$virus.name == virname,]
    }
    # Finally, get rid of rows with NA in them.
    data = data[complete.cases(data),]
    # Check if over filtered.  No point plotting.
    if (nrow(data) == 0) {
        data = NULL
    }
    return(data)
}

# reads in a BPC file which specifies breakpoint coordinates, and filters out data which lies outside of range
# NULL/NA values accepted for chromosome names and ranges, resp.
read.BPC = function(BPC.fn, chrom.A.req, range.A.start, range.A.end, chrom.B.req, range.B.start, range.B.end) {
    data = read.table(BPC.fn, sep='\t', col.names=c("chrom.A", "pos.A", "chrom.B", "pos.B"))

    if (! is.null(chrom.A.req)) data = data[data$chrom.A %in% chrom.A.req,]
    if (! is.null(chrom.B.req)) data = data[data$chrom.B %in% chrom.B.req,]

    if (!is.na(range.A.start) & !is.na(range.A.end)) {
        data = data[data$pos.A >= range.A.start & data$pos.A <= range.A.end,]
    } else if (!is.na(range.A.start)) {
        data = data[data$pos.A >= range.A.start,]
    } else if (!is.na(chrom.range.end)) {
        data = data[data$pos.A <= range.A.end,]
    }

    if (!is.na(range.B.start) & !is.na(range.B.end)) {
        data = data[data$pos.B >= range.B.start & data$pos.B <= range.B.end,]
    } else if (!is.na(range.B.start)) {
        data = data[data$pos.B >= range.B.start,]
    } else if (!is.na(chrom.range.end)) {
        data = data[data$pos.B <= range.B.end,]
    }

    # Get rid of rows with NA in them.
    data = data[complete.cases(data),]

# this is old and should be updated to chrom.A -> chrom.B
#    # Create direction string, listing (human, virus) strand
#    # recall, these are V->H reads
#    data$strand = sprintf("h%sv%s", ifelse(data$chrom.fwd, "+", "-"), ifelse(data$virus.fwd, "+", "-"))

    # Currently unused
    data$chrom.fwd = NULL
    data$virus.fwd = NULL

    return(data)
}


# Note that BreakpointCompare/B_BreakpointPlots.PBP/src/BreakpointPlotter.R uses RP BED files (e.g. PBP_0).  However, we need both virus and human
# PBP positions so use the RP file directly. 
getPindel = function(pindel.fn, is.virus, chrom, chrom.range.start, chrom.range.end) {
    pindel.all = read.table(pindel.fn, col.names=c("chrom1", "start1", "end1", "strand1", "length1", "chrom2", "start2", "end2", "strand2", "length2", 'event.size', "support"), sep="\t", 
                                colClasses = c("character", 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric', 'numeric', 'character', 'numeric', 'numeric', 'character'))
    # gi|310698439|ref|NC_001526.2|   3182    4384    +   1202    14  68683066    68684167    -   0       Support: 45
    # Get rid of reads which aren't either human-virus or virus-human
    pindel.all = pindel.all[xor(grepl('^gi', pindel.all$chrom1), grepl('^gi', pindel.all$chrom2)),]
    # We expect that all chrom1 will be virus and all chrom2 will be human.  Test to make sure this is true.  Other pindel runs may have different parameters
    # and we want to make sure we don't lose data
    if (!all(grepl('^gi', pindel.all$chrom1))) {
        stop("Non-virus discordant reads detected in column 1 of pindel results.")
        # If you get this, see if maybe all virus reads are in chrom2; if not, then need to figure out per row which is human and which is virus
    }
    # calculate 'pos' which is the midpoint between start and end, useful for evaluating whether this section is within region of interest
    pindel.all$pos1 = round(pindel.all$start1 + (pindel.all$end1-pindel.all$start1)/2)
    pindel.all$pos2 = round(pindel.all$start2 + (pindel.all$end2-pindel.all$start2)/2)

    # Filter by chrom
    if (!is.null(chrom)) pindel.all = pindel.all[pindel.all$chrom2 == chrom, ]
    if (!is.na(chrom.range.start)) pindel.all = pindel.all[pindel.all$pos2 > chrom.range.start, ]
    if (!is.na(chrom.range.end)) pindel.all = pindel.all[pindel.all$pos2 <= chrom.range.end, ]
    if (is.virus) {
        pindel = pindel.all[,c("chrom1", "start1", "end1", "pos1")]
    } else {
        pindel = pindel.all[,c("chrom2", "start2", "end2", "pos2")]
    }
    names(pindel) = c("chrom", "start", "end", "pos")

    if (nrow(pindel) == 0) {
        pindel = NULL
    } else {
        pindel$virus = NULL
        pindel=unique(pindel)
    }
    return(pindel)
}

# Return human-human breakpoints as reported in Pindel HHRP.  From BreakpointPlotter.R
read.HHRP = function(HHRP.fn, chrom, chrom.range.start, chrom.range.end) {
    empty.HHRP = data.frame(chrom=character(), pos=integer(), start=integer(), end=integer())
    if (is.null(HHRP.fn)) return(empty.HHRP)
    # based on C_ReadDepth/src/DepthCruncher.R:getPindel
    HHRP.all = read.table(HHRP.fn, col.names=c("chrom1", "start1", "end1", "strand1", "length1", "chrom2", "start2", "end2", "strand2", "length2", 'event.size', "support"), sep="\t",
                               colClasses = c("character", 'numeric', 'numeric', 'character', 'numeric', 'character', 'numeric', 'numeric', 'character', 'numeric', 'numeric', 'character'))
    # gi|310698439|ref|NC_001526.2|   3182    4384    +   1202    14  68683066    68684167    -   0       Support: 45
    # Get rid of reads which connect to virus
    HHRP.all = HHRP.all[!grepl('^gi', HHRP.all$chrom1) & !grepl('^gi', HHRP.all$chrom2),]
    if (nrow(HHRP.all) == 0) {
        return(empty.HHRP)
    } 
    # In the plot we want to use color to indicate paired human-human breakpoints.  To do this, we indicate such pairs with matching
    # values of the variable 'g' (which is also used by CBS for a different purpose).  This then specifies the color of HHRP segments
    HHRP.all$g = rownames(HHRP.all)

    # Create two separate data frames, one for each segment, and combine them.
    HHRP1 = HHRP.all[,c("chrom1", "start1", "end1", "strand1", "g")]
    HHRP2 = HHRP.all[,c("chrom2", "start2", "end2", "strand2", "g")]
    names(HHRP1) = c("chrom", "start", "end", "strand", "g")
    names(HHRP2) = c("chrom", "start", "end", "strand", "g")
    HHRP = rbind(HHRP1, HHRP2)

    # calculate 'pos' which is the midpoint between start and end
    HHRP$pos = round(HHRP$start + (HHRP$end-HHRP$start)/2)

    # finally, remove all entries which lie outside region of interest
    if (!is.null(chrom)) 
        HHRP = HHRP[HHRP$chrom == chrom,]
    if (!is.na(chrom.range.start) && !is.na(chrom.range.end)) {
        keeper = which(HHRP$pos >= chrom.range.start & HHRP$pos <= chrom.range.end)
        HHRP = HHRP[keeper,]
    }

    HHRP=unique(HHRP)

    return(HHRP)
}

# copied from BreakpointPlotter.R
read.rSBP = function(rSBP.fn, is.virus, chrom, chrom.range.start, chrom.range.end) {
    #h.rname h.bp.rpos   v.rname v.bp.rpos
    #14  68683566    gi|310698439|ref|NC_001526.2|   3979
    if (is.null(rSBP.fn))
        return(data.frame(chrom=c(NA), pos=c(NA), chrom2=c(NA), pos2=c(NA)))
    rSBP = read.table(rSBP.fn, sep="\t", header=TRUE, colClasses = c("character", 'numeric', 'character', 'numeric'))
    names(rSBP) = c("chrom", "pos", "chrom2", "pos2")
    if (nrow(rSBP) == 0) return(data.frame(chrom=character(), pos=integer(), chrom2=character(), pos2=integer()))

    if (!is.null(chrom)) 
        rSBP = rSBP[rSBP$chrom == chrom,]
    if (!is.na(chrom.range.start) && !is.na(chrom.range.end)) {
        keeper = which(rSBP$pos >= chrom.range.start & rSBP$pos <= chrom.range.end)
        rSBP = rSBP[keeper,]
    }
    if (is.virus) {
        rSBP$pos = rSBP$pos2
        rSBP$pos2 = NULL
        rSBP$chrom = NULL
    }
    return(rSBP)
}

normalizeCopyNumber = function(read.depth, num.reads, read.length) {
    # Based on CombinedPlot/CombinedPlot.old/orig2013/src/depth_plot.R
    # given read depth data, normalize by average read depth to obtain copy number
    # num_reads is the number of reads in BAM (either total or mapped, your call), typically obtained from .flagstat file
    # read_length is the average read length, also obtained from BAM
    # we normalize the read depth by the average read depth across the genome, defined as 2 * num_reads * read_length / genome_length
    #  - factor 2 is so that copy number comes out to be 2 on average
    num.bp.genome = 3137144693.
    depth.factor = num.reads * read.length / num.bp.genome
#    print(sprintf("Average BAM read depth = %f", depth.factor))
    read.depth$plot.depth = 2. * read.depth$depth / depth.factor
    return(read.depth)
}

# load read depth, filter by chrom and range as necessary, and normalize.  There are three options for normalization:
# - if doLog=TRUE, read.depth$log.depth = log2(depth/<depth>), where depth is read depth and <depth> is the median depth of the filtered dataset
# - if doLog=FALSE:
#   * if num.reads and read.length are specified, norm.depth is read depth divided by average read depth across the genome (as obtained from number of reads and read length)
#   * if num.reads or read.length not specified, depth is returned unchanged.
getReadDepth = function(depth.fn, doLog, chrom=NULL, range.start=NA, range.end=NA, num.reads=NA, read.length=NA) {
    read.depth = read.table(depth.fn, header=FALSE, sep="\t", colClasses=c("character","numeric","numeric"), col.names=c("chrom", "pos", "depth"))

    if (!is.null(chrom)) read.depth = read.depth[read.depth$chrom == chrom, ]
    if (!is.na(range.start)) read.depth = read.depth[read.depth$pos > range.start, ]
    if (!is.na(range.end)) read.depth = read.depth[read.depth$pos <= range.end, ]

    if (doLog) {
        # evaluate log.depth as log2(depth/<depth>)
        read.depth$plot.depth = log2(read.depth$depth/median(read.depth$depth,na.rm=T))
    } else if (!is.na(num.reads) & (!is.na(read.length))) {
        read.depth = normalizeCopyNumber(read.depth, num.reads, read.length)
    } else {
        read.depth$plot.depth=read.depth$depth
    }
    # might be smart to check nrow(read.depth) > 0.  If it is zero, there is nothing to plot.
    return(read.depth)
}

# here we should determine what the y axis label is.
getYaxisLabel = function(doLog, num.reads, read.length) {
    if (doLog) {
        label="Log Depth Ratio"
    } else if (!is.na(num.reads) & (!is.na(read.length))) {
        label="copy number"
    } else {
        label="read depth" 
    }
    return(label)
}

getCBS = function(read.depth, do.log) {
    # Construct CBS data frame which segments read depth into discrete levels
    CBS = calculateCBS(read.depth, do.log)
    names(CBS)[names(CBS)=="loc.start"] = 'start' 
    names(CBS)[names(CBS)=="loc.end"] = 'end'
#      chrom loc.start  loc.end num.mark seg.mean g data.type
#      1    14  68649196 68654404      373  -0.5861       1       CBS     # outdated
    CBS$g = rownames(CBS)  # This becomes the group used for geom_segment, so they are disconnected.
    CBS$data.type = "CBS"
    CBS$num.mark = NULL
    CBS$start = as.numeric(CBS$start)
    CBS$end = as.numeric(CBS$end)
    CBS$pos = round(CBS$start + (CBS$end-CBS$start)/2)
    return(CBS)
}

# TODO: implement read.BPC, analogously to /Users/mwyczalk/Data/TCGA_SARC/ICGC/BreakpointSurveyor/E_Breakpoint/src/BreakpointCruncher.R

# reads in a BPC file which specifies breakpoint coordinates, and filters out data which lies outside of range
# NULL/NA values accepted for chromosome names and ranges, resp.
read.BPC = function(BPC.fn, chrom.A.req, range.A.start, range.A.end, chrom.B.req, range.B.start, range.B.end) {
    data = read.table(BPC.fn, sep='\t', col.names=c("chrom.A", "pos.A", "chrom.B", "pos.B"))

    if (! is.null(chrom.A.req)) data = data[data$chrom.A %in% chrom.A.req,]
    if (! is.null(chrom.B.req)) data = data[data$chrom.B %in% chrom.B.req,]

    if (!is.na(range.A.start) & !is.na(range.A.end)) {
        data = data[data$pos.A >= range.A.start & data$pos.A <= range.A.end,]
    } else if (!is.na(range.A.start)) {
        data = data[data$pos.A >= range.A.start,]
    } else if (!is.na(chrom.range.end)) {
        data = data[data$pos.A <= range.A.end,]
    }

    if (!is.na(range.B.start) & !is.na(range.B.end)) {
        data = data[data$pos.B >= range.B.start & data$pos.B <= range.B.end,]
    } else if (!is.na(range.B.start)) {
        data = data[data$pos.B >= range.B.start,]
    } else if (!is.na(chrom.range.end)) {
        data = data[data$pos.B <= range.B.end,]
    }

    # Get rid of rows with NA in them.
    data = data[complete.cases(data),]

# this is old and should be updated to chrom.A -> chrom.B
#    # Create direction string, listing (human, virus) strand
#    # recall, these are V->H reads
#    data$strand = sprintf("h%sv%s", ifelse(data$chrom.fwd, "+", "-"), ifelse(data$virus.fwd, "+", "-"))

    # Currently unused
    data$chrom.fwd = NULL
    data$virus.fwd = NULL

    return(data)
}

# assume breakpoint.BPC.fn is not NULL
get.discordant.reads = function(breakpoint.BPC.fn, is.virus, chrom, chrom.range.start, chrom.range.end) {
    discordant.reads = read.discordant(breakpoint.BPC.fn)
        # read.BPC = function(BPC.fn, chrom.A.req, range.A.start, range.A.end, chrom.B.req, range.B.start, range.B.end) {
    if (is.virus)
        discordant.reads = discordant.reads[,c("virus.name", "virus.pos")]
    else {
        discordant.reads = discordant.reads[,c("chrom.name", "chrom.pos")]
        if (!is.null(chrom)) 
            discordant.reads = discordant.reads[discordant.reads$chrom.name == chrom,]
        if (!is.na(chrom.range.start) && !is.na(chrom.range.end)) {
            keeper = which(discordant.reads$chrom.pos >= chrom.range.start & discordant.reads$chrom.pos <= chrom.range.end)
            discordant.reads = discordant.reads[keeper,]
        }
    }
    names(discordant.reads) = c("chrom", "pos")
    return(discordant.reads)
}

# legend:
# light gray circles: read depth
# horizontal red lines: CBS
# light gray vertical bars - pindel
# navy blue circles - discordant reads
# multi-colored vertical bars - HHRP
# green dots - rSBP

render.depth = function(read.depth, BPC.df, pindel, CBS, HHRP, rSBP) {
    do.discordant = (!is.null(BPC.df)) && (nrow(BPC.df) > 0)
    do.pindel = !is.null(pindel) && nrow(pindel) > 0
    do.CBS = (!is.null(CBS)) && (nrow(CBS) > 0)
    do.HHRP = !is.null(HHRP) && nrow(HHRP) > 0
    do.rSBP = !is.null(rSBP) && nrow(rSBP) > 0

    # this is to make discordant scatter dots look reasonable.  Questionable whether this is the best way to go tho
    max.y = max(read.depth$plot.depth)
    BPC.df$midpt = max.y/2

    p = ggplot()
    # depth can be one of depth, copy number, or log normalized depth.  read.depth$plot.depth unifies them, and we label the y axes accordingly
    p = p + geom_point(data=read.depth, aes(x=pos, y=plot.depth), alpha=1.0, size=1) 
    #p = p + geom_point(data=data[data$data.type == 'read.depth',], aes(x=pos, y=log.depth), alpha=0.1, size=1) 
    if (do.CBS) 
        p = p + geom_segment(data=CBS, aes(x=start, xend=end, y=seg.mean, yend=seg.mean, group=g), color="#E41A1C", alpha=1) 
        #p = p + geom_segment(data=data[data$data.type == 'CBS',], aes(x=start, xend=end, y=seg.mean, yend=seg.mean, group=g), color="red", alpha=1) 
    if (do.pindel) 
        p = p + geom_rect(data=pindel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="gray50", color="gray50", size=0.25, alpha=0.5)
        #p = p + geom_rect(data=pindel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="gray25", color="gray25", size=0.25, alpha=0.5)
    if (do.discordant)
        p = p + geom_point(data=BPC.df, aes(x=pos, y=midpt), color="#377EB8", size=2.5, alpha=0.25, position = position_jitter(height = max.y/2, width=0))
        #p = p + geom_point(data=data[data$data.type == 'drp',], aes(x=pos, y=0), color="navyblue", size=2.5, alpha=0.25, position = position_jitter(height = 2, width=0))
    if (do.HHRP)
        p = p + geom_rect(data=HHRP, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=g, color=g), size=0.25, alpha=0.25)
    if (do.rSBP)
        p = p + geom_vline(data=rSBP, aes(xintercept=pos, y=0), color="#4DAF4A", alpha=0.5) #, size=4, position = position_jitter(height = 1, width=0))
        #p = p + geom_point(data=rSBP, aes(x=pos, y=0), color="green4", alpha=0.5, size=4, position = position_jitter(height = 1, width=0))

    return(p)
}

annotate.plot = function(p, isVirus, chrom, chrom.range.start, chrom.range.end, ylabel) {
    p = p + theme_bw()
    p = p + theme(legend.position="none")
    p = p + ylab(ylabel)
    if (!isVirus) 
        p = p + xlim(chrom.range.start, chrom.range.end)
    return(p)
}


args = parse_args()

# rule: data will be NULL (and not plotted) if fn not passed, and empty if after filtering no rows.  This may not be yet fully enforced.

#depth.df = get.plot.data(args$depth.fn, args$pindel.fn, args$HHRP.fn, args$rSBP.fn, args$breakpoint.BPC.fn, args$chrom, args$chrom.range.start, args$chrom.range.end, args$virname, args$skip.cbs, args$is.virus)

# pindel.fn, rSBP.fn, HHRP.fn, breakpoint.BPC.fn can all be NULL to skip loading and processing 
# depth.fn is required
# if is.virus=TRUE, extract virus data from pindel, tigra, and discordant reads, otherwise extract from human.
# read.depth, pindel, tigra, and discordant.read.pair data filtered by chrom (this can be e.g. "3" or "gi|310698439|ref|NC_001526.2|"), range.start and range.end

if (!args$is.virus) {
    read.depth = getReadDepth(args$depth.fn, args$plot.log.depth, args$chrom, args$chrom.range.start, args$chrom.range.end, args$num.reads, args$read.length)
} else {
    read.depth = getReadDepth(args$depth.fn, args$plot.log.depth, num.reads=args$num.reads, read.length=args$read.length)
}

if (!is.null(args$breakpoint.BPC.fn)) {
    BPC.df = get.discordant.reads(args$breakpoint.BPC.fn, args$is.virus, args$chrom, args$chrom.range.start, args$chrom.range.end)
} else {
    BPC.df = NULL
} 

if (!is.null(args$pindel.fn)) {
    pindel = getPindel(args$pindel.fn, args$is.virus, args$chrom, args$chrom.range.start, args$chrom.range.end) 
} else {
    pindel = NULL
}

if (!args$skip.cbs) {
# TODO: do this only if log2
    CBS = getCBS(read.depth, args$plot.log.depth)
} else {
    CBS = NULL
}

if (!is.null(args$HHRP.fn)) {
    HHRP = read.HHRP(args$HHRP.fn, args$chrom, args$chrom.range.start, args$chrom.range.end) 
} else {
    HHRP = NULL
}

if (!is.null(args$rSBP.fn)) {
    rSBP = read.rSBP(args$rSBP.fn, args$is.virus, args$chrom, args$chrom.range.start, args$chrom.range.end) 
} else {
    rSBP = NULL
}

depth.ggp = render.depth(read.depth, BPC.df, pindel, CBS, HHRP, rSBP)
ylabel=getYaxisLabel(args$plot.log.depth, args$num.reads, args$read.length) 
depth.ggp = annotate.plot(depth.ggp, args$is.virus, args$chrom, args$chrom.range.start, args$chrom.range.end, ylabel)

if (!args$save.pdf) {
# now save discordant.ggp to file
# http://stat.ethz.ch/R-manual/R-devel/library/base/html/save.html
    cat(paste("Saving to GPP file", args$ggp.fn, "\n"))
    saveRDS(depth.ggp, file=args$ggp.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
# Recall, ggp files can be converted to PDF with ~/bin/ggp2pdf
} else {
    cat(paste("Saving to PDF file", args$ggp.fn, "\n"))
    ggsave(filename=args$ggp.fn, plot=depth.ggp, height=6, width=8, useDingbats=FALSE)
    unlink("Rplots.pdf")
}
