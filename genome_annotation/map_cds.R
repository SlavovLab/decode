
## reconstruct exon lengths to map protein sites to genomic
## positions

## read gff3 exons and 5'UTR information

mam.path <- Sys.getenv("MAMDATA")


## required input
cds.file <- file.path(mam.path, "processedData", "protein_cds_lengths.tsv")
fas.file <- file.path(mam.path, "originalData",
                      "Homo_sapiens.GRCh38.pep.all.fa.gz")
out.file <- file.path(mam.path, "processedData", "protein_cds_structure.dat")
flen <- read.delim(file.path(mam.path, "originalData",
                            "Homo_sapiens.GRCh38.pep.all_lengths.tsv"),
                  header=FALSE, row.names=1)

## load SORTED! cds file
cds <- read.delim(cds.file, header=FALSE, sep=" ")

## list of CDS for each
cdl <- split(cds[,2], cds[,1])

## revert order for minus strand
str <- split(cds[,3], cds[,1])
str <- unlist(lapply(str, unique))
cdl[str=="-"] <- lapply(cdl[str=="-"], rev)

## subtract last CDS position - STOP codon
cdl <- lapply(cdl, function(x) {x[length(x)] <- x[length(x)]-3; x})

## cumsum to get exon positions, last/3 is protein length
cdl <- lapply(cdl, cumsum)

## write out
sink(file=out.file)
for ( i in seq_along(cdl) )
    cat(paste0(names(cdl)[i], "\t", paste(cdl[[i]],collapse=";"), "\n"))
sink()

## TEST

## protein length inferred from CDS reconstruction
mlen <- unlist(lapply(cdl, tail, 1))/3

## protein lengths obtained from FASTA file
rownames(flen) <- sub("\\.[0-9]", "", rownames(flen))
flen <- flen[names(cdl),1]
names(flen) <- names(cdl)

if ( FALSE ) {
    fas <- readFASTA(fas.file, grepID=TRUE)
    names(fas) <- sub("\\.[0-9]", "", names(fas))
    fas <- fas[names(cdl)]
    flen <- lapply(fas, function(x) nchar(x$seq))
}

## get only common
cids <- intersect(names(flen), names(mlen))
mlen <- mlen[cids]
flen <- unlist(flen[cids])
str <- str[cids]

png(file.path(mam.path,"processedData","annotation","protein_cds_qc.png"),
    width=6, height=3, units="in",res=100)
par(mfcol=c(1,2), mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plot(mlen, flen, xlab="protein length inferred from CDS",
     ylab="protein length in fasta file")
abline(a=0, b=1, col=2)
hist(mlen-flen, xlab="length differences", main=NA)
dev.off()

## NOTE: not all CDS sum up to full length protein, off by 1/3 or 2/3
#which(len!=mlen)
## TODO: load transcript fasta and compare to protein sequence

wid <- which.max(abs(mlen-flen))
flen[wid]
mlen[wid]
cdl[[wid]]
if ( exists("fas") )
    fas[wid]
