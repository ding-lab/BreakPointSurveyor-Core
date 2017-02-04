#BPR and BPC file formats

In general, we represent breakpoints as a coordinate given by a pair of
chrom/pos (Breakpoint Coordinates, or BPC).  Alternatively, we may consider a
breakpoint region (BPR), which has a pair of chrom/pos.start/pos.end values.
Breakpoints with precise positions (e.g., discordant pair positions) will be
represented by the former, while regions such as SV Events will be represented
by the latter.

BPC and BPR files may optionaly have a final "attributes" column.  This is a
text field whose interpretation is application-specific. For instance, the BPC
attribute column may code for the color of a point plotted at the given
coordinates.

Each breakpoint coordinate or region is represented just once, with chromA <
chromB (or posA < posB if chromA==chromB)
* Both are represented in a TSV file,
```
       BPC: chromA, posA, chromB, posB, [attribute]
       BPR: chromA, posA.start, posA.end, chromB, posB.start, posB.end, [attribute]
```
* Both datatypes have optional attribute column

BPC and BPR files have no headers.  Lines starting with # are ignored.

#GGP file format

A GGP file is a ggplot object which is saved in binary format.  The utility
`bin/ggp2pdf` converts GGP to PDF for examination

BreakpointSurvey plots are assembled from data in three main steps:

1. Data from a given source is processed ("drawn") to create a GGP file
2. Additional layers can be added to GGP files to incorporate data from additional sources.
3. Multiple GGP files are assembled into one composite BreakpointSurveyor figure and saved as a PDF file.

