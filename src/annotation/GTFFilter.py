#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
#
# Simple script to read GTF file line by line, test if criteria are met, and either print or
# discard line.  Criteria for genes and exons, based on GTF columns:
# Genes: feature = "gene" and source = "protein_coding"
# Exons: feature = "exon" and source = "protein_coding"
#
# In both cases we require the column "source" of the GTF to be protein_coding.  The contents of
# this column (in Ensembl 75 at least) are is inconsistent with GTF documentation.  This is discussed
# here: https://www.biostars.org/p/120306/ Our filtering works for Ensembl 75 GTF but may change in the future.  

import sys
import re

class GTFfeature:
    """Processes one line of GTF file.  Expands attributes column into key/value pairs"""
    # GTF fields described in README.ensembl and http://www.ensembl.org/info/website/upload/gff.html
    def __init__(self, GTFline):
        self.GTF = GTFline
        # TODO: all of these should be inserted into dictionary so can be uniformly accessed.
        self.seqname, self.source, self.feature, start, end, self.score, self.strand, self.frame, attributes = GTFline.split("\t")
        self.start, self.end = int(start), int(end)
        # Split the attributes field by semicolon, then pattern match to extract key, value pair for each.
        attr_list = filter(None, attributes.rstrip().split(";"))
        pattern = re.compile("(\w+?)\s+?\"(.+?)\"")  # matches, 'key "value"'
        # now create key/value dictionary which is saved with this object
        self.attributes = dict(pattern.search(tok).groups() for tok in attr_list)
        # TODO: catch error and print out useful information about,
        #   AttributeError: 'NoneType' object has no attribute 'groups'

    def asGTF(self):
        return self.GTF

    def asBED(self, mergeStrand):
        """Return contents in 6-column BED format: exon, start, end, gene, score, strand
           Start is 0-based
           if mergeStrand, add ":+" or ":-" to gene name 
           """
        gene = self.attributes["gene_name"]
        if mergeStrand: gene += ":"+self.strand
        return "\t".join([self.seqname, str(self.start-1), str(self.end), gene, self.score, self.strand, '\n'])

class GTFexon(GTFfeature):
    def isFeature(self):    # logically, isExon()
        """We consider a line an exon if 1) feature = exon and 2) source = protein_coding"""
        # TODO: this test changes with various releases.  Allow to be specified in a config file.
        #if self.feature == "exon" and self.attributes['gene_biotype'] == "protein_coding": return True
        if self.feature == "exon" and self.source == "protein_coding": return True
        return False

class GTFgene(GTFfeature):
    def isFeature(self):  # logically, isGene()
        """We consider a line a gene if 1) feature = gene and 2) source = protein_coding"""
        #if self.feature == "gene" and self.attributes['gene_biotype'] == "protein_coding": return True
        if self.feature == "gene" and self.source == "protein_coding": return True
        return False

def filterGTF(f, o, GTFclass, asBED, mergeStrand):
    """Go through GTF file line by line and evaluate - using filter method of class specified by the calling function - whether to 
       retain this line or not.  If line retained, output it in either GTF or BED format"""
    for l in f:
        if l.startswith("#"): continue 
        gtf = GTFclass(l)
        if gtf.isFeature():   # true for either genes or exons, depending on what class we were passed.
            if asBED: o.write(gtf.asBED(mergeStrand))
            else: o.write(gtf.asGTF())

def main():
    from optparse import OptionParser
    usage_text = """usage: %prog [options] feature ...
        Read an ENSEMBL GTF file and retain select features, writing output in GTF or BED format.
        'feature' is either 'gene' or 'exon'.  This corresponds to GTF 'feature' column of that name and 'gene_biotype' attribute of 'protein_coding'
        """

    parser = OptionParser(usage_text, version="$Revision: 1.0 $")
    parser.add_option("-v", dest="verbose", action="store_true", help="verbose output")
    parser.add_option("-f", dest="in_fn", default="stdin", help="GTF filename")
    parser.add_option("-o", dest="out_fn", default="stdout", help="Output filename")
    parser.add_option("-b", dest="as_bed", action="store_true", help="Output in BED format")
    parser.add_option("-s", dest="mergeStrand", action="store_true", help="Add :-delimited strand information to gene name")

    (options, params) = parser.parse_args()
    if options.verbose:
        sys.stderr.write(str(options)+'\n')
        sys.stderr.write(str(params)+'\n')

    if len(params) != 1:
        parser.error("Please pass feature name (gene or exon).")

    feature = params[0]
    if feature not in ["gene", "exon"]:
        parser.error("Feature must be 'gene' or 'exon'.")

    if options.in_fn == "stdin":
        f = sys.stdin 
    else:
        f = open(options.in_fn, 'r')

    if options.out_fn == "stdout":
        o = sys.stdout
    else:
        o = open(outfn, "w")

    if feature == "gene":
        filterClass = GTFgene
    else:
        filterClass = GTFexon
    filterGTF(f, o, filterClass, options.as_bed, options.mergeStrand)

    o.close()

if __name__ == '__main__':
    main()
