# Contig 

Creating SBP files ("SAM BreakPoint") from bwa-mem-aligned reads.  These contain
base-pair resolution information about breakpoints detected between chimeric segment.

1. [SAMReader.py](SAMReader.py): takes SAM file, finds breakpoints for every line.
        (note, contigs span two lines).  Writes pSBP file, which has the columns,
'''    "query_name", "ref_name", "bp_pos", "end_pos", "ref_pos", "is_reverse", "is_left"'''
    (ref_pos is refrerence position of nearest mapped position to breakpoint)
2. [BreakPointParser.R](BreakPointParser.R): Creates SBP file by merging a pair of pSBP lines based on query_name
    to create the SBP.  Segments Sa and Sb correspond to left and right end of chimeric
    segment, resp. Multiple breakpoints per contig result in multiple (N-1) lines in BreakpointParser
3. [SBPprocessor.R](SBPprocessor.R): creates simple BPC (aka rSBP = reduced SBP) file with chromA/chromB
    breakpoint positions (1-index format).  A qSBP file is optionally
    generated which gives information about paired breakpoints (more than one
    breakpoint per contig)

Tigra-sv sometimes creates pathologically long contig names 
which makes samtools (and pysam) choke.  [qname_convert.py](src/contig/qname_convert.py)
shortens QNAMEs in all SAM files to a hex string using an MD5 hash.

