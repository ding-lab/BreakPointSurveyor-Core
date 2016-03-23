#!/usr/bin/python
# Matthew A. Wyczalkowski, mwyczalk@genome.wustl.edu
#
# Simple script to read GTF file line by line, test if criteria are met, and either print or
# discard line.  Criteria for genes and exons are based on GTF documentation and guesswork.
#
# The first version Ensembl this was developed for is Ensembl 75, which had an inconsistency between
# the documentation and implementation (see https://www.biostars.org/p/120306/).  
# We define the following ad hoc criteria to extract exon and gene position, name, and strand information:
#   Ensembl 75:
#       Genes: feature = "gene" and source = "protein_coding"
#       Exons: feature = "exon" and source = "protein_coding"
# By Ensembl 84 the formatting of the GTF file had changed, and the following criteria are being used
# to extract gene and exon positions:
#   Ensembl 84:
#       Genes: feature = "gene"
#       Exons: feature = "exon" and source = "ensembl"
#       (the source criterion seems to limit to  
# The version of Ensembl being used can be passed with the -e flag.  Discussion in the link above suggets that 
# versions 77 and later have this issue fixed and can be parsed like ensembl 84.
#
# TODO: it would be useful to allow user to define criteria for gene/exon exclusion by passing arguments to this program.

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

class GTFexon75(GTFfeature):
    def isFeature(self):    # logically, isExon()
        """Evaluate feature from ensembl 75 GTF file.  Format of this version seems inconsistent with documentation, 
        and our approach is ad hoc.
        We consider a line an exon if 1) feature = exon and 2) source = protein_coding"""
        if self.feature == "exon" and self.source == "protein_coding": return True
        return False

class GTFgene75(GTFfeature):
    def isFeature(self):  # logically, isGene()
        """Evaluate feature from ensembl 75 GTF file.  Format of this version seems inconsistent with documentation, 
        and our approach is ad hoc.
        We consider a line a gene if 1) feature = gene and 2) source = protein_coding"""
        if self.feature == "gene" and self.source == "protein_coding": return True
        return False

class GTFexon84(GTFfeature):
    def isFeature(self):    
        """Evaluate feature from ensembl 84 GTF file.  Other releases may also work.
        We consider a line an exon if 1) feature = exon and 2) source = ensembl.
        This is ad hoc. 
        """
        if self.feature == "exon" and self.source == "ensembl": return True
        return False

class GTFgene84(GTFfeature):
    def isFeature(self):  # logically, isGene()
        """Evaluate feature from ensembl 84 GTF file.  Other releases may also work.
        We consider a line a gene if feature = gene.
        This is ad hoc."""
        if self.feature == "gene": return True
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
    parser.add_option("-e", dest="ensembl_version", default="84", help="Define Ensembl version.  Supported values = '75', '84'.")

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

    if options.ensembl_version not in ['75', '84']:
        print "Unsupported Ensembl version "+options.ensembl_version
        print "Assuming format is Ensembl 84"
        options.ensembl_version = '84'

    if feature == "gene":
        filterClass = GTFgene75 if options.ensembl_version == '75' else GTFgene84
    else:
        filterClass = GTFexon75 if options.ensembl_version == '75' else GTFexon84
    filterGTF(f, o, filterClass, options.as_bed, options.mergeStrand)

    o.close()

if __name__ == '__main__':
    main()
