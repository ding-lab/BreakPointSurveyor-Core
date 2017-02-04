# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
#
# Usage: Rscript PvalBubblePlotter.R [-v] [-t title] [-l] [-p pval] [-P logP.max] [-F] [-h height] [-w width] [-x x-label text]
#                -a int.start -b int.end -f flank -s case.sample.id gene.pval.dat out.pdf
# Visualize gene dysregulation in vicinity of integration event
# 
# -v verbose output
# -a, -b: start, end of integration site 
# -f specify "flank" for plot, the distance on either side of integration event in the figure.
# -s: the case sample ID (e.g., TCGA-BA-4077-01)
# -t title: title of plot
# -x axis.text.x: X axis label [NULL]
# -l: hide legend
# -p: for each bubble plot a guide mark at given p-value
# -P logP.max: maximum value of -log10(P); larger values will be cropped to this size.  Default 10
# -F: plot FDR instead of P-value.  Requires FDR column to be present in data.  All p-value-related parameters (e.g., logP.max) will refer to FDR instead.
# -h, -w: plot height and width in inches [default 2, 8, resp.]

library("plyr")
suppressMessages(library("ggplot2"))
library("RColorBrewer")
library("grid")

library(scales)

options("width"=180) # useful for debugging
options(scipen=3)  # no scientific notation

#print("Will quit on warning.")  # do this only for debugging warning messages
#options(warn=2)


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
    # named arguments.  Note that Rscript will not process '-g' correctly - ugh.
    verbose = get_bool_arg(args, "-v")
    int.start = as.numeric(get_val_arg(args, "-a", NA))
    int.end = as.numeric(get_val_arg(args, "-b", NA))
    flank = as.numeric(get_val_arg(args, "-f", NA))
    case.sample.id = get_val_arg(args, "-s", NULL)
    pval.guide = as.numeric(get_val_arg(args, "-p", NA))
    logP.max = as.numeric(get_val_arg(args, "-P", 10))
    title = get_val_arg(args, "-t", NULL)
    axis.title.x = get_val_arg(args, "-x", NULL)
    hide.legend = get_bool_arg(args, "-l")
    plot.FDR = get_bool_arg(args, "-F")
    plot.height = as.numeric(get_val_arg(args, "-h", 2))
    plot.width = as.numeric(get_val_arg(args, "-w", 8))
    
    out.fn = args[length(args)];             args = args[-length(args)]
    gene.pval.fn = args[length(args)];             args = args[-length(args)]


    val = list( 'verbose' = verbose, 'int.start' = int.start, 'int.end' = int.end, 'flank' = flank, 'case.sample.id' = case.sample.id,
                'title' = title, 'out.fn' = out.fn, 'gene.pval.fn' = gene.pval.fn, 'hide.legend'=hide.legend, 'pval.guide'=pval.guide, 'logP.max'=logP.max,
                'plot.FDR'=plot.FDR, 'plot.height'=plot.height, 'plot.width'=plot.width, 'axis.title.x'=axis.title.x)
    if (val$verbose) { print(val) }
    return (val)
}

# This is based closely on file below; see there for annotation code and other useful tidbits.
# /Users/mwyczalk/Data/Virus/Virus_2013.9a/RSEM-Exon/RPKM-Scatter/src/RPKM_scatter_plotter.R

# read in p-value data and calculate logP.  If plot.FDR, logP is the corrected (FDR) p-value.
get.gene.pval = function(data.fn, plot.FDR) {
    data = read.table(data.fn, header=TRUE)
#                     integration.event          gene     stream strand chrom    start      end p.value up.regulated
#1  TCGA-BA-4077-01B-01D-2268-08.chr14A         ACTN1   upstream      -    14 69340859 69446157 0.04621         TRUE

    # calculate the midpoint of gene
    data$gene.mid = data$start + (data$end - data$start) / 2

    # define arrow.start, arrow.end based on strand of gene
    forward = data$strand == "+"
    data$arrow.start = ifelse(forward, data$start, data$end)
    data$arrow.end = ifelse(forward, data$end, data$start)

    if (plot.FDR) {
        data$logP = -log10(data$FDR)
    } else {
        data$logP = -log10(data$p.value)
    }
    
    return(data)
}

# Create bubble plot visualizing P value for genes in vicinity of breakpoint.  pval guides serve to indicate pvalue relative to a reference
# bubble size indicates -ln(P.value) and whether up or down regulated 
# text labels indicate genes, gene size and direction indicated by arrows.  Note that genes with integration events within them are split,
# and left, right exon groups are indicated as separate genes.
plot.Pval = function(Pval.data, int.start, int.end, range.start, range.end, title, logP.guide, logP.max, hide.legend=FALSE) {
    # logP.crop becomes the minimum of logP and logP.max, so that limit() does not delete data
    Pval.data$logP.max = logP.max
    Pval.data[is.infinite(Pval.data$logP),"logP"] = logP.max
    Pval.data$logP.crop = apply(Pval.data[,c("logP", "logP.max")], 1, min) 

    Pval.data$logP.guide = logP.guide

    Pval.data$annotation.pos = 0.7
    Pval.data$arrow.pos = 0.75
    Pval.data$hjust = 1

    reverse.gene = which(Pval.data$strand == "-")
    if (length(reverse.gene) > 0) {
        Pval.data[reverse.gene,"annotation.pos"] = -Pval.data[reverse.gene,"annotation.pos"]   ###
        Pval.data[reverse.gene,"arrow.pos"] = -Pval.data[reverse.gene,"arrow.pos"] 
        Pval.data[reverse.gene,"hjust"] = 0
    }

    Pval.data$regulated = ifelse(Pval.data$up.regulated, "Upregulated", "Downregulated")
    color.hex = c("#E41A1C", "#377EB8") #Set 1 colors, starting with red:
    color.names = c("Upregulated", "Downregulated")
    names(color.hex) = color.names
    gene.color.scale = scale_color_manual(name="regulated", values=color.hex)#, guide=TRUE)

    # Gene position is plotted at midpoint
    p = ggplot(data=Pval.data) 

    # Switching data here so are not overplotting the rectangle repeatedly
    p = p + geom_rect(data=data.frame(foo=0), alpha=0.5, xmin=int.start, xmax=int.end, ymin=-1, ymax=5, fill="#4DAF4A", color="#4DAF4A", size=0.25)  

    myarrow = arrow(length=unit(0.1,"cm"), type="closed")
    p = p + geom_segment(mapping=aes(x=arrow.start, xend=arrow.end, y=arrow.pos, yend=arrow.pos), arrow=myarrow, color='gray50', size=0.25, alpha=0.25, show.legend=FALSE)

    p = p + geom_point(aes(x=gene.mid, y=0, color=regulated, size=logP.crop), alpha = 0.5)  # V2

    if (!is.na(logP.guide))
        p = p + geom_point(aes(x=gene.mid, y=0, size=logP.guide), color="gray60", shape=1, alpha=0.25)

    p = p + gene.color.scale
    breaks.scale = -log10(c(0.05, 5e-5, 5e-8))
    breaks.labels = c("P=0.05", "P=5e-5", "P=5e-8")
    p = p + scale_size_area(limit=c(0,logP.max), breaks=breaks.scale, labels=breaks.labels, max_size=20)

    p = p + geom_text(aes(x=gene.mid, y=annotation.pos, label=gene, hjust=hjust), angle=90, size=2.5, alpha = 0.75, color="gray50")


    p = p + theme_bw() + theme(panel.background = element_blank(), panel.border = element_blank())

    #p = p + scale_color_brewer(palette="Set1")

    xmin = min(Pval.data$start, range.start)
    xmax = max(Pval.data$end, range.end)

    p = p + scale_x_continuous(labels=comma, limits=c(xmin,xmax)) + ylim(-1,1)  # comma requires library(scales)

    if (is.null(args$axis.title.x)) {
        p = p + theme(axis.title.x = element_blank())
    } else {
        p = p + xlab(args$axis.title.x) 
        p = p + theme(axis.title = element_text(size=8, color="gray50"))
    }
    p = p + theme(legend.title=element_blank(), axis.ticks = element_blank())
    p = p + theme(axis.text.y = element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())

    p = p + ggtitle(title) + theme(plot.title = element_text(size=8, color="gray50"))
    p = p + theme(axis.text = element_text(size=8, color="gray50"))
    if (hide.legend) {
        p = p + theme(legend.position="none", axis.title.y=element_blank())
    }
    return(p)
}



args = parse_args()

Pval.data = get.gene.pval(args$gene.pval.fn, args$plot.FDR)

range.start = max(0, args$int.start - args$flank)
range.end = args$int.end + args$flank

if (!is.na(args$pval.guide)) {
    logP.guide = -log10(args$pval.guide)
} else {
    logP.guide = NA
}

Pval.ggp = plot.Pval(Pval.data, args$int.start, args$int.end, range.start, range.end, args$title, logP.guide, args$logP.max, args$hide.legend)

cat(sprintf("    Saved to %s\n", args$out.fn))
ggsave(filename=args$out.fn, Pval.ggp, height=args$plot.height, width=args$plot.width, useDingbats=FALSE)
unlink("Rplots.pdf")

# other stuff
