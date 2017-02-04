# BreakPointSurveyor-Core
There are three layers of BreakPoint Surveyor project:

* BPS Core: core analysis and plotting, typically in R or Python
* BPS Workflow: Project- and locale-specific workflows. Mostly as BASH scripts
* BPS Data: BPS-generated secondary data, graphical objects, and plots

BreakPoint Surveyor is currently released with reference implementations of the Workflow and Data layers for the TCGA_Virus project:
* [BreakPointSurveyor-Workflow.TCGA_Virus](https://github.com/ding-lab/BreakPointSurveyor-Workflow.TCGA_Virus)
* [BreakPointSurveyor-Data.TCGA_Virus](https://github.com/ding-lab/BreakPointSurveyor-Data.TCGA_Virus)

Multi-panel figures are generated in three steps: 

1. The data processing normalizes data into standard formats. For instance, breakpoint
predictions from different SV callers are normalized into a [BPC](doc/FileFormat.md) file format,
while read depth and gene annotation are converted to Depth and BED formats, respectively.  
2. Each dataset is rendered as an image panel saved as a binary "GGP" object.
Additional layers, for instance predictions from different SV callers, may be added to an existing
GGP object in subsequent processing steps.
3. Finally, multiple GGP objects are assembled, aligned to common axes,
and saved to a PDF format to form a composite figure.

Scripts
-------
### src/analysis
*Utilities for RPKM expression analysis, read depth, Pindel output processing*

**[ExonExpressionAnalyzer.R](src/analysis/ExonExpressionAnalyzer.R)**
Evaluates relative gene expression based on RPKM data from case and control.
Calculate and write to stdout p-value associated with gene expression in vicinity of integration event.
**TODO** describe algorithm

**[ExonPicker.R](src/analysis/ExonPicker.R)**
Select exons from genes upstream and downstream of integration event and write BED file describing these.

**[Pindel_RP.Reader.R](src/analysis/Pindel_RP.Reader.R)**
Create Breakpoint Region file ([BPR](doc/FileFormat.md)) based on output of Pindel RP module.  

**[RPKM_Joiner.R](src/analysis/RPKM_Joiner.R)**
Process multiple RPKM files and combine column-wise into one data file.  

**[TigraCTXMaker.R](src/analysis/TigraCTXMaker.R)**
Create a breakdancer-style CTX file from either Pindel's RP or [BPR](doc/FileFormat.md) data to be used as Tigra-SV input

**[depthFilter.py](src/analysis/depthFilter.py)**
Read BAM file and evaluate read depth in a segment. Output is subsampled to give data size,
optimized for performance.

**[vafFilter.py](src/analysis/vafFilter.py)**
Parse VAF as output by Pindel

### src/annotation
*Ad hoc scripts for processing Ensembl gene/exon names and regions.*

**[ChromRenamer.py](src/annotation/ChromRenamer.py)**
Translate chromosome names in BED file between two standards using a database.
Used for normalizing feature names, as discussed [here](https://www.biostars.org/p/138011/)

**[GTFFilter.py](src/annotation/GTFFilter.py)**
Simple script to read GTF file line by line, test if criteria are met, and either print or
discard line.  Used to extract gene and exon domains.

**[TLAExamine.R] (src/annotation/TLAExamine.R)**
Tool for examining GTF and VCF files.  Expands a column of key/value pairs into multiple columns,
for examining in e.g. spreadsheet.

### src/contig
*Utilities related to contig creation with Tigra-SV*

Contig alignment improves breakpoint predictions by assembling a consensus
sequence (contig) from reads spanning a breakpoint, then re-aligning the contig
to the human+virus reference.  Contigs are created using Tigra-SV based on
breakpoint predictions from Pindel and aligned with BWA mem (v. 0.7.10) to the
human+virus reference described above. The resulting SAM file is then processed
to yield a series of files (described below) used for plotting and analysis,
with support for multiple breakpoints per contig.

* pSBP: One entry per alignment (line) of SAM file, describing alignment of
  one segment to reference
* SBP: Incorporates pair of pSBP segments, Sa and Sb, describing one
  breakpoint.  Multiple (N) breakpoints in a contig will result in (N-1) SBP
  lines
* BPC (aka rSBP): Retains only chromosome names and breakpoint position
* qSBP: Gives information about paired breakpoints (>1 breakpoints per contig).
  For N breakpoints per contig, N-1 qSBP entries are generated.  The BPC file
  is used in the breakpoint panel of the BPS structure plot to illustrate
  contig predictions.

Creating SBP files ("SAM BreakPoint") from bwa-mem-aligned reads.  These contain
base-pair resolution information about breakpoints detected between chimeric segment.

1. [SAMReader.py](src/contig/SAMReader.py): takes SAM file, finds breakpoints for every line.
        (note, contigs span two lines).  Writes pSBP file, which has the columns,
'''    "query_name", "ref_name", "bp_pos", "end_pos", "ref_pos", "is_reverse", "is_left"'''
    (ref_pos is refrerence position of nearest mapped position to breakpoint)
2. [BreakPointParser.R](src/contig/BreakPointParser.R): Creates SBP file by merging a pair of pSBP lines based on query_name
    to create the SBP.  Segments Sa and Sb correspond to left and right end of chimeric
    segment, resp. Multiple breakpoints per contig result in multiple (N-1) lines in BreakpointParser
3. [SBPprocessor.R](src/contig/SBPprocessor.R): creates simple BPC (aka rSBP = reduced SBP) file with human, 
    virus breakpoint positions (1-index format).  A qSBP file is optionally
    generated which gives information about paired breakpoints (more than one
    breakpoint per contig)

Tigra-sv sometimes creates pathologically long contig names 
which makes samtools (and pysam) choke.  [qname_convert.py](src/contig/qname_convert.py)
shortens QNAMEs in all SAM files to a hex string using an MD5 hash.

### src/plot
*Heart of BPS figure generation*

Each dataset is rendered as an image panel using the ggplot() function and
saved as a binary "GGP" object with saveRDS().
GGP objects can be visualized using the `ggp2pdf` utility.  Additional layers,
for instance predictions from different SV callers, may be added to an existing
GGP object in a subsequent processing step with data from a different BPC (or
BPR) file. 

**[AnnotationDrawer.R](src/plot/AnnotationDrawer.R)**
Create gene annotation GGP files with optional exon definitions for each gene.

**[BreakpointDrawer.R](src/plot/BreakpointDrawer.R)**
Common BreakpointSurveyor plotting utilities.

**[BreakpointSurveyAssembler.R](src/plot/BreakpointSurveyAssembler.R)**
Create or append various features to breakpoint coordinate GGP file.  Chrom A coordinates are plotted
on X axis, B on Y.

**[DepthDrawer.R](src/plot/DepthDrawer.R)**
Plot read depth (or related quantites) over a genomic region and add annotation to this plot.

**[DepthUtil.R](src/plot/DepthUtil.R)**
Common read depth utilities.

**[HistogramDrawer.R](src/plot/HistogramDrawer.R)**
Create a histogram of read depth (or estimated copy number) for chrom A and B

**[PvalBubblePlotter.R](src/plot/PvalBubblePlotter.R)**
Visualize gene dysregulation in vicinity of integration event

**[ZoomGGP.R](src/plot/ZoomGGP.R)**
Utility to change plot limits of a GGP file and save as PDF

**[ggp2pdf](src/plot/ggp2pdf)**
Convert GGP (ggplot binary) file into pdf

### src/util
*Common and ad hoc scripts*

**[BPS_Util.R](src/util/BPS_Util.R)**
Common BreakpointSurveyor utilities.

**[PlotListMaker.py](src/util/PlotListMaker.py)**
Create a Breakpoint Surveyor PlotList file from Breakpoint Coordinate (BPC) or Breakpoint Region (BPR) data

**[PlotListParser.R](src/util/PlotListParser.R)**
Given barcode, chrom, and chrom position, return PlotList name which contains this position.

**[makeBreakpointRegions.py](src/util/makeBreakpointRegions.py)**
Cluster Breakpoints into Breakpoint regions.

**[processVCF.py](src/util/processVCF.py)**
Read VCF file and write coordinates of features in various formats
