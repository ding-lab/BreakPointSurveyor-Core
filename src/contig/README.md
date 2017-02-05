# Contig 

## Background and Workflow
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

## Utilities
Creating SBP files ("SAM BreakPoint") from bwa-mem-aligned reads.  These contain
base-pair resolution information about breakpoints detected between chimeric segment.

1. [SAMReader.py](SAMReader.py): takes SAM file, finds breakpoints for every line.
        (note, contigs span two lines).  Writes pSBP file, which has the columns,
'''    "query_name", "ref_name", "bp_pos", "end_pos", "ref_pos", "is_reverse", "is_left"'''
    (ref_pos is refrerence position of nearest mapped position to breakpoint)
2. [BreakPointParser.R](BreakPointParser.R): Creates SBP file by merging a pair of pSBP lines based on query_name
    to create the SBP.  Segments Sa and Sb correspond to left and right end of chimeric
    segment, resp. Multiple breakpoints per contig result in multiple (N-1) lines in BreakpointParser
3. [SBPprocessor.R](SBPprocessor.R): creates simple BPC (aka rSBP = reduced SBP) file with human, 
    virus breakpoint positions (1-index format).  A qSBP file is optionally
    generated which gives information about paired breakpoints (more than one
    breakpoint per contig)

Tigra-sv sometimes creates pathologically long contig names 
which makes samtools (and pysam) choke.  [qname_convert.py](src/contig/qname_convert.py)
pdf(file=args$out.fn, width=args$width, height=args$height, useDingbats=FALSE)
shortens QNAMEs in all SAM files to a hex string using an MD5 hash.
