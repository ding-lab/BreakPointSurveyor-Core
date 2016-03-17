# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript AnnotationRenderer.R [-v] [-P] [-V] [-A range] [-e exons.bed] [-D dodge] [-B] genes.bed out.ggp
#
# Create gene annotation GGP files with optional exon definitions for each gene.  
#
# Input arguments:
#
# * genes.bed is a genes BED file.  Genes outside of the range given by (chrom, range.start, range.end) will be discarded.
# * out.ggp is filename of GGP output  If -P defined, .pdf appended to filename if necessary

# Optional arguments:
# -A: Define range, specified as "C" or "C:M-N", where C is chromosome name
#         and M,N are genomic start, stop positions.  Range is used for two distinct purposes:
#           1) Filter data before plotting
#           2) Specify plot region
#         Filtering data may be prevented with -F.  
# -D annotation dodge parameter.  Annotations will fall into this many lines.
# -e exons.bed.  BED file which indicates regions of exons.  Similar to genes.bed.
# -B: Annotate for chrom B.  This rotates text and reverses X scale.
# -P: Output as PDF file instead of GGP.  This is primarily for convenience and debugging.

options("width"=180) # useful for debugging
library("bitops")
library("plyr")
suppressPackageStartupMessages(library("ggplot2"))
#require(gtable)

source_relative = function(source.fn) {
    # http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
    initial.options = commandArgs(trailingOnly = FALSE)
    file.arg.name = "--file="
    script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename = dirname(script.name)
    other.name = paste(sep="/", script.basename, source.fn)
#    print(paste("Sourcing",other.name,"from",script.name))
    source(other.name)
}
source_relative("BPS_Util.R")

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    pdf.out = get_bool_arg(args, "-P")

    range.A = parse.range.str(get_val_arg(args, "-A", "all"))  # accessible as range.chr, start=range.pos[1], start=range.pos[2]
    exons.bed.fn = get_val_arg(args, "-e", NULL)
    dodge = get_val_arg(args, "-D", "4")
    is.B = get_bool_arg(args, "-B")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.ggp = args[length(args)];             args = args[-length(args)]
    genes.bed.fn = args[length(args)];      args = args[-length(args)]

    val = list( 'range.pos'=range.A$range.pos, 'range.chr'=range.A$range.chr,  'verbose' = verbose, 
            'exons.bed.fn' = exons.bed.fn, 'dodge' = as.numeric(dodge), 'pdf.out'=pdf.out,
            'out.ggp' = out.ggp, 'genes.bed.fn' = genes.bed.fn, 'is.B'=is.B)
    if (val$verbose) { print(val) }

    return (val)
}


# Annotation data frames have the following columns:
# xmax, xmin, ymax, ymin - determines position of annotation box
# label - text annotation.  May be blank.
# labelx, labely - label text positions, relevant only when label not blank.
# label_group - determines color of annotation.  Same value as label when label exists, but should not be blank.
#
# these are then plotted directly with make.annotation.ggp
# we make sure labels are unique (no genes listed twice) 

# Creates annotation data frame which describes genes to annotate
# If no annotation, zero-row data frame with a subset of columns returned
get.chrom.annotation.df = function(genes, range.pos, dodge, height=1) {
    annotation = data.frame(xmin=genes$start, xmax=genes$end, label=genes$name, label_group=genes$name)
    # if no annotation, return data frame now.
    if (nrow(annotation) == 0) return(annotation)
    # in some cases a gene is listed twice (e.g. KRBOX1).  To keep gene names unique, collapse all duplicates
    # (identified by label and label_group) and keep the superset of positions
    annotation = ddply(annotation, c("label", "label_group"), summarise, xmin = min(xmin), xmax = max(xmax))

    # we define range.min, range.max as the ranges of the drawn genes, considering that they are cut off at range start, end
    annotation$range.start = range.pos[1]
    annotation$range.end = range.pos[2]
    annotation$range.min = apply(annotation[,c("xmin", "range.start")],1,max)  # FYI, pmax/pmin is cleaner way to do this
    annotation$range.max = apply(annotation[,c("xmax", "range.end")],1,min)

    # label is positioned in middle of drawn region of gene
    annotation$labelx=annotation$range.min+(annotation$range.max-annotation$range.min)/2   
    # now get rid of unnecessary columns
    annotation$range.start = NULL
    annotation$range.end = NULL
    annotation$range.min = NULL
    annotation$range.max = NULL

    # For legibility in the plot, y position of each gene is given by ypos = i %% dodge, where i is the rank of position along X axis
    # this gives genes a "stairstep" layout
    annotation$ymin = (rank(annotation$labelx)-1) %% dodge  
    annotation$ymax = annotation$ymin + height
    annotation$labely = annotation$ymax + 0.3
    return(annotation)
}

get.gene.exon.annotation.df = function(genes, exons, range.pos, dodge) { 
    annotation.df = get.chrom.annotation.df(genes, range.pos, dodge, 0.1)  # last number is "thickness" of intron section 
    if (nrow(exons) != 0 & nrow(annotation.df) != 0) {
        # define position of exons so they line up with gene y positions
        exon_annotation = data.frame(xmin=exons$start, xmax=exons$end, label="", label_group=exons$name)

        # y position of exon is based on y position of gene
        # In certain cases we find there are exons with no corresponding genes.  Rather than crashing, we warn the user of this
        # and go on.  Deal this with all.x in merge, and detect it by finding exon_annotation with unknown ymin, which comes from 
        # gene annotation.
        exon_annotation = merge(exon_annotation, annotation.df[,c("label_group", "ymin", "labelx", "labely")], by="label_group", sort=FALSE, all.x=TRUE) 
        if (any(is.na(exon_annotation$ymin))) {
            is.na.ymin = which(is.na(exon_annotation$ymin))
            unknown.genes = unique(exon_annotation[is.na.ymin, "label_group"])
            cat(paste("Warning: these genes have exon data but are not on gene list: ", paste(unknown.genes), "\n"))
            exon_annotation[is.na.ymin, "ymin"] = 0
        }
        exon_annotation$ymin = exon_annotation$ymin - 0.25   
        exon_annotation$ymax = exon_annotation$ymin + 0.6   # this number corresponds to height of exon

        annotation.df = rbind(annotation.df, exon_annotation)
    }
    return (annotation.df)
}

# Create GGP object from annotation data frame. 
# is.B indicats that this is a panel which will be rotated when assembled in final plot.  This implies,
#   * text angle is rotated 90 degrees
#   * X axis is reversed
make.annotation.ggp = function(annotation, is.B) {    
    # text.angle is the angle of the annotation text next to the genes
    text.angle = if (is.B) -90 else 0

    ggp = ggplot() 
    ggp = ggp + geom_rect(data=annotation, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=label_group, color=label_group), alpha=0.1) 

    vjust=-0.2
    ggp = ggp + geom_text(data=annotation, aes(x=labelx, y=labely, label=label, color=label_group), size=3, vjust=vjust, angle=text.angle) 
    ggp = ggp + theme(legend.position="none")
    ggp = ggp + theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), 
         plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
         panel.background = element_blank())
    ggp = ggp + xlab(NULL) + ylab(NULL) 
    ggp = ggp + ylim(min(annotation$ymin),max(annotation$ymax)+1)  # make enough room for the labels.
    if (is.B)
        ggp = ggp + scale_x_reverse()  # make x position increase in the downward direction for B panel only
    return(ggp)
}

args = parse_args()

genes = read.BED(args$genes.bed.fn)
genes = filter.BED(genes, args$range.chr, args$range.pos)

exons = read.BED(args$exons.bed.fn)
exons = filter.BED(exons, args$range.chr, args$range.pos)

# all genes and exons within range will be drawn
annotation.df = get.gene.exon.annotation.df(genes, exons, args$range.pos, args$dodge)

# If there is no annotation we simply don't write a file and let assembly deal with this
if (nrow(annotation.df)==0) {
    cat(paste("No annotation after filtering.  Not saving",args$out.ggp,"\n"))
    q()
}

ggp = make.annotation.ggp(annotation.df, args$is.B)

# now save ggp to file
write.GGP(ggp, args$out.ggp, args$pdf.out)
