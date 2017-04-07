# BreakPointSurveyor-Core
*Core utilities implementing BreakPointSurveyor workflow.*

## Overview

BreakPointSurveyor (BPS) is a set of core libraries (this project) and 
[workflows](https://github.com/ding-lab/BreakPointSurveyor) which, with optional external tools,
evaluate genomic sequence data to discover, analyze, and provide a visual summary of
interchromosomal breakpoint events.

The [BreakPointSurveyor project](https://github.com/ding-lab/BreakPointSurveyor) provides three reference workflows, each implemented as a separate [git
branch](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell).
These workflows (and the links to view them) are:

* **TCGA_Virus** ([`master` branch](https://github.com/ding-lab/BreakPointSurveyor)):  Comprehensive workflow and data for one
  TCGA virus-positive sample
  ([TCGA-BA-4077-01B-01D-2268-08](https://gdc-portal.nci.nih.gov/legacy-archive/files/6533e56c-b5b8-4c85-862b-a5526c5c2e0a))
  which has been aligned to a custom reference
* **1000SV** ([`1000SV` branch](https://github.com/ding-lab/BreakPointSurveyor/tree/1000SV)): Analysis of discordant reads on publicly available human sample 
* **Synthetic** ([`Synthetic` branch](https://github.com/ding-lab/BreakPointSurveyor/tree/Synthetic)): Creation and analysis of a dataset containing an inter-chromosomal breakpoint 

**Citation** *In prep*


## Getting Started

See [BreakPointSurveyor documentaton](https://github.com/ding-lab/BreakPointSurveyor) and 
the [installation instructions](https://github.com/ding-lab/BreakPointSurveyor/blob/master/INSTALL.md).


## Documentation
### Architecture
There are three layers of BreakPointSurveyor (BPS) project:

* BPS Core: core analysis and plotting, typically in R or Python
* BPS Workflow: Project- and locale-specific workflows. Mostly as BASH scripts
* BPS Data: BPS-generated secondary data, graphical objects, and plots

The [BreakPointSurveyor](https://github.com/ding-lab/BreakPointSurveyor) provides
three example workflows and their data.  This project ([BreakPointSurveyor-Core](https://github.com/ding-lab/BreakPointSurveyor-Core))
provides the Core layer.  It is typically distributed as a submodule of the [BreakPointSurveyor](https://github.com/ding-lab/BreakPointSurveyor)
project and does not need to be installed separately.

### Visualization
Multi-panel figures are generated in three steps: 

1. The data processing normalizes data into standard formats. For instance, breakpoint
predictions from different SV callers are normalized into a [BPC](https://github.com/ding-lab/BreakPointSurveyor/blob/master/Development.md)) file format,
while read depth and gene annotation are converted to Depth and BED formats, respectively.  
2. Each dataset is rendered as an image panel saved as a binary "GGP" object.
Additional layers, for instance predictions from different SV callers, may be added to an existing
GGP object in subsequent processing steps ([see details](https://github.com/ding-lab/BreakPointSurveyor/blob/master/Development.md)).
3. Finally, multiple GGP objects are assembled, aligned to common axes,
and saved to a PDF format to form a composite figure.


## BPS Utilities

BPS Core consists of a number of utilities which are used by Workflow scripts to process and visualize data.  They are described
below, ordered by directory structure.

### src/analysis
*Utilities for RPKM expression analysis, read depth, Pindel output processing*

* **[ExonExpressionAnalyzer.R](src/analysis/ExonExpressionAnalyzer.R)**
Evaluates relative gene expression based on RPKM data from case and control.
Calculate and write to stdout p-value associated with gene expression in vicinity of integration event.
[Algorithm details](https://github.com/ding-lab/BreakPointSurveyor/blob/master/L_Expression/AlgorithmDetails.md).

* **[ExonPicker.R](src/analysis/ExonPicker.R)**
Select exons from genes upstream and downstream of integration event and write BED file describing these.

* **[Pindel_RP.Reader.R](src/analysis/Pindel_RP.Reader.R)**
Create Breakpoint Region file ([BPR](https://github.com/ding-lab/BreakPointSurveyor/blob/master/Development.md))) based on output of Pindel RP module.  

* **[RPKM_Joiner.R](src/analysis/RPKM_Joiner.R)**
Process multiple RPKM files and combine column-wise into one data file.  

* **[TigraCTXMaker.R](src/analysis/TigraCTXMaker.R)**
Create a breakdancer-style CTX file from either Pindel's RP or [BPR](https://github.com/ding-lab/BreakPointSurveyor/blob/master/Development.md)) data to be used as Tigra-SV input

* **[depthFilter.py](src/analysis/depthFilter.py)**
Read BAM file and evaluate read depth in a segment. Output is subsampled to give data size,
optimized for performance.

* **[vafFilter.py](src/analysis/vafFilter.py)**
Parse VAF as output by Pindel

### src/annotation
*Ad hoc scripts for processing Ensembl gene/exon names and regions.*

* **[ChromRenamer.py](src/annotation/ChromRenamer.py)**
Translate chromosome names in BED file between two standards using a database.
Used for normalizing feature names, as discussed [here](https://www.biostars.org/p/138011/)

* **[GTFFilter.py](src/annotation/GTFFilter.py)**
Simple script to read GTF file line by line, test if criteria are met, and either print or
discard line.  Used to extract gene and exon domains.

* **[TLAExamine.R] (src/annotation/TLAExamine.R)**
Tool for examining GTF and VCF files.  Expands a column of key/value pairs into multiple columns,
for examining in e.g. spreadsheet.

### src/contig
*Utilities related to contig creation with Tigra-SV*

Contig alignment improves breakpoint predictions by assembling a consensus
sequence (contig) from reads spanning a breakpoint, then re-aligning the contig
to the human+virus reference.  Contigs are created using Tigra-SV.

Read more about [contig workflow here](src/contig/README.md).

### src/plot
*BPS figure rendering and assembly*

Each dataset is rendered as an image panel using the ggplot() function and
saved as a binary "GGP" object with saveRDS().
GGP objects can be visualized using the `ggp2pdf` utility.  Additional layers,
for instance predictions from different SV callers, may be added to an existing
GGP object in a subsequent processing step with data from a different BPC (or
BPR) file. 

* **[AnnotationDrawer.R](src/plot/AnnotationDrawer.R)**
Create gene annotation GGP files with optional exon definitions for each gene.

* **[BreakpointDrawer.R](src/plot/BreakpointDrawer.R)**
Common BreakpointSurveyor plotting utilities.

* **[BreakpointSurveyAssembler.R](src/plot/BreakpointSurveyAssembler.R)**
Create or append various features to breakpoint coordinate GGP file.  Chrom A coordinates are plotted
on X axis, B on Y.

* **[DepthDrawer.R](src/plot/DepthDrawer.R)**
Plot read depth (or related quantites) over a genomic region and add annotation to this plot.

* **[DepthUtil.R](src/plot/DepthUtil.R)**
Common read depth utilities.

* **[HistogramDrawer.R](src/plot/HistogramDrawer.R)**
Create a histogram of read depth (or estimated copy number) for chrom A and B

* **[PvalBubblePlotter.R](src/plot/PvalBubblePlotter.R)**
Visualize gene dysregulation in vicinity of integration event

* **[ZoomGGP.R](src/plot/ZoomGGP.R)**
Utility to change plot limits of a GGP file and save as PDF

* **[ggp2pdf](src/plot/ggp2pdf)**
Convert GGP (ggplot binary) file into pdf

### src/util
*Common and ad hoc scripts*

* **[BPS_Util.R](src/util/BPS_Util.R)**
Common BreakpointSurveyor utilities.

* **[PlotListMaker.py](src/util/PlotListMaker.py)**
Create a Breakpoint Surveyor PlotList file from Breakpoint Coordinate (BPC) or Breakpoint Region (BPR) data

* **[PlotListParser.R](src/util/PlotListParser.R)**
Given barcode, chrom, and chrom position, return PlotList name which contains this position.

* **[makeBreakpointRegions.py](src/util/makeBreakpointRegions.py)**
Cluster Breakpoints into Breakpoint regions.

* **[processVCF.py](src/util/processVCF.py)**
Read VCF file and write coordinates of features in various formats

## Authors
Matthew A. Wyczalkowski, m.wyczalkowski@wustl.edu

## License
This software is licensed under the GNU General Public License v3.0

## Acknowledgements

This work was supported by the National Cancer Institute [R01CA178383 and
R01CA180006 to Li Ding, R01CA172652 to Ken Chen]; and National Human Genome Research
Institute [U01HG006517 to Li Ding]. 

