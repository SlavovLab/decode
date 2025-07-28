
## extract coding sequences from transcript fasta via UTR coordinates


library(segmenTools)
options(stringsAsFactors=FALSE)


#' retrieve codon from a sequence object by position
getCodon <- function(x, pos) {
    sq <- strsplit(x$seq,"")[[1]]
    paste0(sq[pos + 1:3], collapse="")
}

## read gff3 exons and 5'UTR information

mam.path <- Sys.getenv("MAMDATA")
out.path <- file.path(mam.path, "processedData", "annotation")
out.file <- file.path(mam.path, "processedData", "transcript_utr_lengths.tsv")
out.fasta <- file.path(mam.path, "processedData", "coding.fa")

## INPUT FILES
gff.file <- file.path(mam.path,"originalData",
                      "Homo_sapiens.GRCh38.110_utr.gff3.gz")
fasta.file <- file.path(mam.path, "originalData",
                        "Homo_sapiens.GRCh38.cdna.all.fa.gz")
tpmap.file <- file.path(mam.path, "originalData",
                        "protein_transcript_map.tsv")
feature.file <- file.path(mam.path, "features_GRCh38.110.tsv")


## parse gff3 genome annotation file
gff <- gff2tab(gff.file)

## keep only required columns
gff <- gff[,c("feature","start","end","strand","Parent")]

## remove transcript tag
gff$Parent <- sub("^transcript:","", gff$Parent)

## split by transcript!
models <- split(gff, f=gff$Parent)

## sort each gene transcript model
models <- lapply(models, function(x) x[order(x$start,x$end),])

## convert coordinates to local!
models <- lapply(models, function(x) {
    x[,c("start","end")] <-
        x[,c("start","end")] - min(x[,c("start","end")]) +1
    x
})


## CHECK: all gene models are on the same strand!
strands <- lapply(models, function(x) unique(x$strand))
table(unlist(lapply(strands, length)))

## CALCULATE UTR LENGTHs
## NOTE: MULTIPLE UTR if splicing within UTR 
utr5 <- lapply(models, function(x) {
    x <- x[x$feature=="five_prime_UTR",,drop=FALSE]
    sum(x$end-x$start+1) # sums up multiple
})
## same for 3'UTR
utr3 <-  lapply(models, function(x) {
    x <- x[x$feature=="three_prime_UTR",,drop=FALSE]
    sum(x$end-x$start+1)
})

utr3 <- unlist(utr3)
utr5 <- unlist(utr5)

utr <- data.frame(transcript=names(utr5),
                  start=utr5, end=utr3[names(utr5)])

## check mutual availability
## TODO: add 0 if one UTR but not the other is available
if ( sum(!names(utr5)%in%names(utr3))>0 )
    stop("not all 5'UTR have a 3'UTR")
if ( sum(!names(utr3)%in%names(utr5))>0 )
    stop("not all 3'UTR have a 5'UTR")



## LOAD TRANSCRIPT FASTA AND CUT
fas <- readFASTA(fasta.file)
## name sequences by transcript ID
desc <- lapply(fas, function(x) {unlist(strsplit(x$desc, " "))})
ids <- unlist(lapply(desc, function(x) sub("\\..*","",x[1])))
names(fas) <- ids

## TODO: understand missing!
## ENSP00000481152 - ENST00000619729
## ENST00000439082
## ENSP00000354876 - ENST00000361739 
if ( interactive() ) {
    tid="ENST00000361739"
    tid%in%names(fas)    # present: transcript contains
    tid%in%rownames(utr) # not present: no UTR
}

## FILTER: only keep transcripts and UTRs for which all info is available
##cids <- intersect(rownames(utr), names(fas))
##cat(paste("cutting", length(cids), "fasta sequences\n"))
##fas <- fas[rownames(utr)]

## FILTER: only do transcripts that code for a protein
tpmap <- read.delim(tpmap.file, header=FALSE)
cids <- names(fas)%in%tpmap[,1]
cat(paste("writing", sum(cids), "coding sequences to", out.fasta,"\n"))
fas <- fas[cids]


## add 0 UTR if no UTR is available
missutr <- names(fas)[!names(fas)%in%utr[,1]]
if ( length(missutr)>0 ) {
    cat(paste("adding", length(missutr),
              "UTR of length 0 to coding sequences\n"))
    mutr <- cbind.data.frame(transcript=missutr, start=0, end=0)
    utr <- rbind(utr, mutr)

}


if ( !interactive() ) {
    ## write out UTR lengths
    write.table(utr, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
}

if ( !interactive() ) {

    ## write out FASTA
    sink(file=file(out.fasta, open = "wt"))
    for ( i in seq_along(fas) ) {
        id <- names(fas)[i]
        ##id <- utr[i,"transcript"]
        len <- nchar(fas[[id]]$seq)
        start <- 1
        end <- len
        if ( id %in% rownames(utr) ) {
            start <- start + utr[id,"start"]
            end <-   len   - utr[id,"end"]
        }
        cutseq <- substr(fas[[id]]$seq, start, end)
        cat(paste0(">", id, ", coding region, ",
                   start,"-",end, "\n", cutseq, "\n"))
    }
    sink()
}


### SANITY CHECK: LOAD CODING FASTA produced above,
### and do stats over STOP and START codons.
tfas <- readFASTA(out.fasta, grepID=TRUE)
cids <- names(tfas)

## FILTER: only take MANE-tagged!
features <- read.delim(feature.file)
cids <- cids[cids%in%features$MANE]

tfas <- tfas[cids]


fig.name <- file.path(out.path,"utr_mapping_sanity_2.png")
## subset for quicker interactive plot
if ( interactive() & length(cids)>2e4 ) {
    set.seed(1)
    ridx <- sample(ridx, 10000)
    fig.name <- file.path(out.path,"utr_mapping_sanity_test_2.png")
} 

seqlen <- unlist(lapply(tfas, function(x) nchar(x$seq)))
endpos <- seqlen -3
names(endpos) <- names(tfas)
ends   <- sapply(cids, function(id) getCodon(tfas[[id]], endpos[id]))
starts <- sapply(cids, function(id) getCodon(tfas[[id]], 0))

## summarize
N <- length(starts)
startst <- sort(table(starts), decreasing=TRUE)
endst <- sort(table(ends), decreasing=TRUE)

## plot sanity check
png(fig.name, res=300, width=5, height=5, units="in")
par(mfcol=c(2,1), mai=c(.5,1,.15,.1), mgp=c(2.5,.3,0), tcl=-.25, xaxs="i")
barplot(100*startst/N,las=2, log="", ylab="% of start codons")
legend("top", paste(N, "transcripts"))
barplot(100*endst/N,las=2, log="", ylab="% of stop codons")
dev.off()



