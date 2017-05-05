# Matthew Wyczalkowski
# m.wyczalkowski@wustl.edu
# The Genome Institute
#
# Common BreakpointSurveyor plotting utilities.

# save ggp object to either GGP binary representation or PDF file.
# if writing to PDF (pdf.out=TRUE), append ".pdf" to filename out.fn if it does not already have that extension
write.GGP = function(ggp, out.fn, pdf.out=FALSE) {
    # http://stat.ethz.ch/R-manual/R-devel/library/base/html/save.html
    if (!pdf.out) {
        cat(paste("Saving to GPP file", out.fn, "\n"))
        saveRDS(ggp, file=out.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/
    } else {

        ext = file_ext(out.fn)
        if (ext != "pdf")
           out.fn = paste(sep=".", out.fn, "pdf") 

        cat(paste("Saving to PDF file", out.fn, "\n"))

        ggsave(plot=ggp, filename=out.fn, useDingbats=FALSE)
        unlink("Rplots.pdf")
    }
}
