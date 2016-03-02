# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript AnnotationCruncher.R [-v] [-V] [-a range.start] [-b range.end] [-c chrom] [-e exons.bed] [-D dodge] [-B] genes.bed annotation.ggp
#
# Create gene annotation GGP files with optional exon definitions for each gene.  
#
# Input arguments:
#
# genes.bed is a genes BED file.  Genes outside of the range given by (chrom, range.start, range.end) will be discarded.
# annotation.gpp is output file in GGP format

# Mandatory:
# -a range.start, -b range.end, -c chrom: chromosome name and start, end positions of region of interest.  1-based.  Mandatory.
# Optional arguments:
# -D annotation dodge parameter.  Annotations will fall into this many lines.
# -e exons.bed.  Bed file which indicates regions of exons.  Similar to genes.bed.
# -B: rotate text -90 degrees.  Useful for annotating chrom B
#
# v2.0 - virus annotation deprecated.
# v1.1 - fixed annotation bugs based on /Users/mwyczalk/Data/Virus/Virus_2013.9a/RSEM-Exon/RPKM-Scatter/src/RPKM_scatter_plotter.R 


options("width"=180) # useful for debugging
library("bitops")
library("plyr")
suppressPackageStartupMessages(library("ggplot2"))
#require(gtable)

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
    range.start = get_val_arg(args, "-a", NA)  # default should be NA for things converted to numeric
    range.end = get_val_arg(args, "-b", NA)
    chrom = get_val_arg(args, "-c", "all")
    exons.bed.fn = get_val_arg(args, "-e", NULL)
    dodge = get_val_arg(args, "-D", "4")
    label.flip = get_bool_arg(args, "-B")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    ggp.fn = args[length(args)];             args = args[-length(args)]
    genes.bed.fn = args[length(args)];      args = args[-length(args)]

    val = list( 'verbose' = verbose, 'range.start' = as.numeric(range.start),
            'range.end' = as.numeric(range.end), 'chrom' = chrom, 'exons.bed.fn' = exons.bed.fn, 'dodge' = as.numeric(dodge),
            'ggp.fn' = ggp.fn, 'genes.bed.fn' = genes.bed.fn, 'label.flip'=label.flip)
    if (val$verbose) { print(val) }

    return (val)
}

# be able to read bed files of arbitrary number of columns.  Columns named as described here:
#   http://bedtools.readthedocs.org/en/latest/content/general-usage.html
read.bed = function(bed.fn) {
    #data = read.csv(bed.fn, header = FALSE, sep = "\t", col.names=c("chrom", "start", "end", "name"), row.names=NULL)
    # this complains if file has zero lines
    data = read.csv(bed.fn, header = FALSE, sep = "\t", row.names=NULL)
    names(data)[1:3] = c("chrom", "start", "end")
    if (length(names(data)) > 3) names(data)[4] = "name" 
    if (length(names(data)) > 4) names(data)[5] = "score" 
    if (length(names(data)) > 5) names(data)[6] = "strand" 
    return(data)
}

# read bed file with gene or exon definitions
read.filtered.bed = function(bed.fn, chrom, range.start, range.end) {
    data = read.bed(bed.fn)
    # discard rows (exons or genes) which do not overlap with range of interest.
    # Do this by excluding rows whose 1) start is greater than range end or 2) end is less than range start
    data = data[data$chrom == chrom,]
    data = data[data$start < range.end,]
    data = data[data$end > range.start,]
    return(data)
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
get.chrom.annotation.df = function(genes, range.start, range.end, dodge, height=1) {
    annotation = data.frame(xmin=genes$start, xmax=genes$end, label=genes$name, label_group=genes$name)
    # in some cases a gene is listed twice (e.g. KRBOX1).  To keep gene names unique, collapse all duplicates
    # (identified by label and label_group) and keep the superset of positions
    annotation = ddply(annotation, c("label", "label_group"), summarise, xmin = min(xmin), xmax = max(xmax))

    # we define range.min, range.max as the ranges of the drawn genes, considering that they are cut off at range.start, range.end
    annotation$range.start = range.start
    annotation$range.end = range.end
    annotation$range.min = apply(annotation[,c("xmin", "range.start")],1,max)
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

get.gene.exon.annotation.df = function(genes, exons, range.start, range.end, dodge) { 
    annotation.df = get.chrom.annotation.df(genes, range.start, range.end, dodge, 0.1)  # last number is "thickness" of intron section 
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
# text.angle is the angle of the annotation text next to the genes
make.annotation.ggp = function(annotation, text.angle=0) {    
    p = ggplot() 

    # note that cropped annotations are not implemented
    p = p + geom_rect(data=annotation, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=label_group, color=label_group), alpha=0.1) 

    vjust=-0.2
    p = p + geom_text(data=annotation, aes(x=labelx, y=labely, label=label, color=label_group), size=3, vjust=vjust, angle=text.angle) 
    p = p + theme(legend.position="none")
    p = p + theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), 
         plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
         panel.background = element_blank())
    p = p + xlab(NULL) + ylab(NULL) 
    p = p + ylim(min(annotation$ymin),max(annotation$ymax)+1)  # make enough room for the labels.
    return(p)
}

args = parse_args()

exons = read.filtered.bed(args$exons.bed.fn, args$chrom, args$range.start, args$range.end)
genes = read.filtered.bed(args$genes.bed.fn, args$chrom, args$range.start, args$range.end)

# all genes and exons within range will be drawn
annotation.df = get.gene.exon.annotation.df(genes, exons, args$range.start, args$range.end, args$dodge)

# If there is no annotation we simply don't write a file and let assembly deal with this
if (nrow(annotation.df)==0) {
    #annotation.ggp = rectGrob(gp = gpar(col = "white"))  # can't do this in make.a.ggp because coord_flip dies on a rectGrob
    cat(paste("No annotation after filtering.  Not saving",args$ggp.fn,"Downstream workflow will ignore annotation.\n"))
    q()
}

text.angle = if (args$label.flip) -90 else 0
annotation.ggp = make.annotation.ggp(annotation.df, text.angle)

# now save annotation.ggp to file
# http://stat.ethz.ch/R-manual/R-devel/library/base/html/save.html
cat(paste("Saving to GPP file", args$ggp.fn, "\n"))
saveRDS(annotation.ggp, file=args$ggp.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/

