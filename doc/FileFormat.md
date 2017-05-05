# BPR and BPC file formats

We represent breakpoints as a coordinate given by a pair of
chromosome names and positions (breakpoint coordinates, or BPC).  

Similarly, a breakpoint region (BPR), has a pair of chromosome, start, and end values.

BPC and BPR files may optionaly have a final "attributes" column.  This is a
text field whose interpretation is application-specific. For instance, the BPC
attribute column may code for the color of a point plotted at the given
coordinates.

Each breakpoint coordinate or region is represented just once, with chromA <
chromB (or posA < posB if chromA==chromB), in a tab-separated (TSV) format:
```
       BPC: chromA  posA  chromB  posB  [attribute]
       BPR: chromA  posA.start  posA.end  chromB  posB.start  posB.end  [attribute]
```

BPC and BPR files have no headers.  Lines starting with # are ignored.

#GGP file format

A GGP file is a ggplot object (i.e., output of `ggplot()` in
[GGPlot](http://ggplot2.org/) package), saved in binary format.  

Use the utility `ggp2pdf` to convert GGP to PDF:

```bps-core/src/plot/ggp2pdf test.ggp test.pdf``` 

