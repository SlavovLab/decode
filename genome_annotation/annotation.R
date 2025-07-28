library(segmenTools)
options(stringsAsFactors=FALSE)

## util
imap <- function(...) invisible(Map(...))

mam.path <- Sys.getenv("MAMDATA")

## INPUT FILES
gff.file <- file.path(mam.path,"originalData",
                      "Homo_sapiens.GRCh38.110_genes.gff3.gz")

mrna.file <- file.path(mam.path,"originalData",
                       "Homo_sapiens.GRCh38.110_transcripts.gff3.gz")
prot.file <- file.path(mam.path,"originalData",
                       "Homo_sapiens.GRCh38.110_proteins.gff3.gz")

## Transcript <-> Protein mapping
tpmap.file <- file.path(mam.path,"originalData", "protein_transcript_map.tsv")

## ensembl<>refseq mapping
refseq.file <- file.path(mam.path,"originalData",
                         "ensembl_refseq_20240528.tsv.gz")

## gene name synonyms, via https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
## downloaded on 20240712
synonym.file <- file.path(mam.path, "originalData", "gene_synonyms.tsv")

## GO annotation
go.file <- file.path(mam.path,"originalData",
                     "ensembl_gene_goslim_20231204.tsv.gz")
gof.file <- file.path(mam.path,"originalData",
                      "ensembl_gene_go_20231204.tsv.gz")
up.file <- file.path(mam.path,"originalData",
                     "ensembl_gene_uniprot_20231204.tsv.gz")



## OUTPUT FILE
feature.file <- file.path(mam.path, "features_GRCh38.110.tsv")

fig.path <- file.path(mam.path,"processedData","annotation")
dir.create(fig.path, showWarnings=FALSE)
                      

## MAIN GENE TABLE

## parse gff3 genome annotation file
gff <- gff2tab(gff.file)
## required columns from gff3
gff.cols <- c(ID="gene_id", name="Name", type="biotype",
              chr="seqname", start="start", end="end", strand="strand",
              description="description")

genes <- gff[,gff.cols]
colnames(genes) <- names(gff.cols)
rownames(genes) <- genes$ID

cat(paste("loaded", nrow(genes), "genes\n"))

## get transcripts
mrna <- gff2tab(mrna.file)

## any non-transcripts?
grep("transcript:", mrna$ID, invert=TRUE, value=TRUE)

mrna$ID <- sub("transcript:","",mrna$ID)

## get canonical and MANE tags
table(grep("MANE_Select",mrna$tag,value=T))
table(grep("canonical",mrna$tag,value=T))
mrna$canonical <- rep("",nrow(mrna))
mrna$canonical[grep("canonical",mrna$tag)]<- mrna$ID[grep("canonical",mrna$tag)]
mrna$MANE <- rep("",nrow(mrna))
mrna$MANE[grep("MANE_Select",mrna$tag)] <- mrna$ID[grep("MANE_Select",mrna$tag)]

## TRANSCRIPT start/end sites
tssc <- mrna$start
tssc[mrna$strand=="-"] <- mrna$end[mrna$strand=="-"]

tesc <- mrna$end
tesc[mrna$strand=="-"] <- mrna$start[mrna$strand=="-"]
## check correct orientation
table(mrna$strand[tssc<tesc]);table(mrna$strand[tssc>tesc])

## split by parent transcript
tssl <- split(tssc, f=sub("^gene:","",mrna$Parent))
tesl <- split(tesc, f=sub("^gene:","",mrna$Parent))

## GENE start/end sites
starts <- genes$start
starts[genes$strand=="-"] <- genes$end[genes$strand=="-"]
names(starts) <- genes$ID

ends <- genes$end
ends[genes$strand=="-"] <- genes$start[genes$strand=="-"]
names(ends) <- genes$ID



## QC check coordinates start>end for minus strand
table(starts < genes$end, genes$strand)
if ( sum(!names(tssl)%in%names(starts))>0 )
    stop("missing start sites\n")
if ( sum(!names(tesl)%in%names(ends))>0 )
    stop("missing end sites\n")

## get those for which we have start sites
starts <- starts[names(tssl)]
ends <- ends[names(tesl)]
##starts[genes$strand=="-"] <- - starts[genes$strand=="-"]

## RELATIVE TSS: transcript starts  - gene start
## -> should be positive for upstream TSS!
tssl <- Map(`-`, tssl, starts)
## revert for negative strand values
strnd <- genes[names(tssl),"strand"]
tssl[strnd=="-"] <- lapply(tssl[strnd=="-"], function(x) -x)

## RELATIVE TES: transcript end - gene end
## -> should be negative for upstream TES!
tesl <- Map(`-`, tesl, ends)
## revert for negative strand values
strnd <- genes[names(tesl),"strand"]
tesl[strnd=="-"] <- lapply(tesl[strnd=="-"], function(x) -x)

## NOTE: TSS are >=0 and TES are <=0
## --> LONGEST TRANSCRIPT APPEARS TO DEFINE GENE!!


## TRANSCRIPTS by Parent -> gene
mrnl <- split(mrna$transcript_id, f=sub("^gene:","",mrna$Parent))


## QC: duplicated transcript per gene?
table(unlist(lapply(mrnl, function(x) sum(duplicated(x)))))

## get proteins for matched transcripts
prot <- gff2tab(prot.file)

## check: all CDS with same protein ID from unique transcripts?
cat(paste(length(unique(prot$ID)), "unique proteins\n"))
cat(paste(length(unique(prot$Parent)),
          "unique parent transcripts of proteins\n"))

## .. and reduce to unique proteins
prot <- prot[!duplicated(prot$ID),]

## split by transcript parent
prot <- unlist(split(prot$protein_id, f=sub("^transcript:","",prot$Parent)))

## NOTE: only 1 protein per transcript!
table(unlist(lapply(prot,length)))



cat(paste(sum(names(prot)%in%mrna$transcript_id),
          "proteins WITH a transcript parent\n"))
cat(paste(sum(!names(prot)%in%mrna$transcript_id),
          "proteins W/O a transcript parent\n"))

## attach proteins to coding transcript parents 
mrna$proteins <- prot[match(mrna$transcript_id, names(prot))]

## replace MANE and canonical transcripts by their proteins
## NOTE: should only be one!
##mrna$canonical <- prot[mrna$canonical]
##mrna$MANE <- prot[mrna$MANE]
canl <- split(mrna$canonical, f=sub("^gene:","",mrna$Parent))
manl <- split(mrna$MANE, f=sub("^gene:","",mrna$Parent))

## get proteins for the main transcripts
tpmap <- read.delim(tpmap.file, header=FALSE, row.names=1)

manll <- lapply(manl, function(x) x[x!=""])
mapl <- manll
mapl[lengths(manll)>0] <-  tpmap[unlist(manll),1] 


## MAP PROTEINS TO GENES
prol <- lapply(split(mrna$proteins, f=sub("^gene:","",mrna$Parent)),
               function(x) x[!is.na(x)])

cat(paste("transcripts found for", sum(genes$ID%in%names(mrnl)), "genes\n"))
cat(paste("proteins found for", sum(genes$ID%in%names(prol)), "genes\n"))


## GENE NAME SYNONYMS
syns <- read.delim(synonym.file)
synl <- lapply(split(syns[,2], syns[,1]), unique)
synl <- lapply(synl, function(x) unique(unlist(strsplit(x, "\\|"))))
synl <- lapply(synl, function(x) x[!x%in%c("","-")])
synl <- unlist(lapply(synl, paste, collapse=";"))

### ANNOTATIONS 

## GOslim
go <- read.delim(go.file)
## filter empty
go <- go[go[,1]!="",]
## collect and WRITE OUT GO TERMS to be used for later analyses
terms <- go[!duplicated(go[,3]),3:4]
write.table(terms, file=file.path(mam.path,"processedData","goslim.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

## GO list by gene
gol <- split(go[,3], f=go[,1])
gol <- lapply(gol, unique)

## full GO
gof <- read.delim(gof.file)
## filter empty
gof <- gof[gof$GO.term.accession!="",]
## GO list by gene
gofl <- split(gof[,3], f=gof[,1])
gofl <- lapply(gofl, unique)

## UNIPROT IDs
up <- read.delim(up.file)

## SWISSPROT
spl <- split(up$UniProtKB.Swiss.Prot.ID, f=up[,1])
spl <- lapply(spl, unique)
spl <- lapply(spl, function(x) x[x!=""])

## all uniprot
upl <- split(up$UniProtKB.Gene.Name.ID, f=up[,1])
upl <- lapply(upl, unique)
upl <- lapply(upl, function(x) x[x!=""])

## all uniprot isoforms
ipl <- split(up$UniProtKB.isoform.ID, f=up[,1])
ipl <- lapply(ipl, unique)
ipl <- lapply(ipl, function(x) x[x!=""])


## merge all into main genes table

collaps <- function(x) {
    rm <- x%in%c("NULL","NA","") | is.na(x)
    paste0(unique(x[!rm]), collapse=";")
}

## match and unlist

canonical <- unlist(lapply(canl[match(genes$ID, names(mrnl))], collaps))
mane <- unlist(lapply(manl[match(genes$ID, names(mrnl))], collaps))
mane.protein <- unlist(lapply(mapl[match(genes$ID, names(mrnl))], collaps))
transcripts <- unlist(lapply(mrnl[match(genes$ID, names(mrnl))], collaps))
TSS <- unlist(lapply(tssl[match(genes$ID, names(tssl))],collaps))
TES <- unlist(lapply(tesl[match(genes$ID, names(tesl))],collaps))
proteins <- unlist(lapply(prol[match(genes$ID, names(prol))], collaps))
swissprot <- unlist(lapply(spl[match(genes$ID, names(spl))], collaps))
uniprot <-   unlist(lapply(upl[match(genes$ID, names(upl))], collaps))
isoprot <-   unlist(lapply(ipl[match(genes$ID, names(ipl))], collaps))
go <- unlist(lapply(gol[match(genes$ID, names(gol))], collaps))
gof <- unlist(lapply(gofl[match(genes$ID, names(gofl))], collaps))

## match by name
synl <- synl[match(genes$name, names(synl))]

agenes <- cbind(genes,
                synonyms=synl,
                canonical=canonical,
                MANE=mane,
                MANE.protein=mane.protein,
                swissprot=swissprot,
                transcripts=transcripts,
                TSS=TSS,
                TES=TES,
                proteins=proteins,
                uniprot=uniprot,
                uniprot.isoforms=isoprot,
                GOslim=go, GO=gof)


### WRITE OUTPUT FILE                                              
write.table(agenes,
            file=feature.file,
            sep="\t", row.names=FALSE, quote=FALSE, na="")


### STATISTICS and PLOTS
## statistics table and figures
## transcripts, proteins, UTRs per gene etc.


W <- 5
H <- 3.5
ftyp <- "png"
res <- 200

tgenes <- read.delim(feature.file)

## restrict to protein coding
tgenes <- tgenes[tgenes$type=="protein_coding",]

## transcription start sites

tssl <- lapply(strsplit(tgenes$TSS, ";"),as.numeric)
tesl <- lapply(strsplit(tgenes$TES, ";"),as.numeric)

## get exceptionally long distances
tsslng <- unlist(lapply(tssl, max))
teslng <- unlist(lapply(tesl, min))
tss.long <- tsslng >  1500e3
tes.long <- teslng < -1500e3


plotdev(file.path(fig.path, "TES_distance"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unlist(tesl)/1e3, xlab="TES distance from `gene` end / kb",
     breaks=50, main="transcript end sites")
text(teslng[tes.long]/1e3, 1000,
     labels=tgenes$name[tes.long], srt=90, pos=4, cex=.7)
dev.off()
plotdev(file.path(fig.path, "TSS_distance"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unlist(tssl)/1e3, xlab="TSS distance from `gene` start / kb",
     breaks=50, main="transcript start sites")
text(tsslng[tss.long]/1e3, 1000,
     labels=tgenes$name[tss.long], srt=90, pos=4, cex=.7)
dev.off()
## zoom

plotdev(file.path(fig.path, "TES_distance_zoom"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unlist(tesl)/1e3, xlab="TES distance from `gene` end / kb",
     xlim=c(-100,0), breaks=1000, main="transcript end sites")
dev.off()
plotdev(file.path(fig.path, "TSS_distance_zoom"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unlist(tssl)/1e3, xlab="TSS distance from `gene` start / kb",
     xlim=c(0,100), breaks=1000, main="transcript start sites")
dev.off()


## transcripts per gene
mrna <- unlist(lapply(strsplit(tgenes$transcripts, ";"),length))
plotdev(file.path(fig.path, "transcripts_per_gene"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(mrna, breaks=0:max(mrna),
     xlab="transcripts per `gene`",
     main=paste(sum(mrna), "ensembl transcripts"))
legend("topright", c(paste(" 1 transcript /gene:",sum(mrna==1)),
                     paste(" 2 transcripts/gene:",sum(mrna==2)),
                     paste(">2 transcripts/gene", sum(mrna>2))), bty="n")
text(mrna[mrna>50], 100, labels=tgenes$name[mrna>50], srt=90, pos=4, cex=.7)
dev.off()

## ENSP proteins per gene
unil <- unlist(lapply(strsplit(tgenes$proteins, ";"),length))
plotdev(file.path(fig.path, "proteins_per_gene"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unil, breaks=0:max(unil),
     xlab="ensembl protein IDs per `gene`",
     main=paste(sum(unil), "ensembl proteins"))
legend("topright", c(paste(" 1 protein /gene:",sum(unil==1)),
                     paste(" 2 proteins/gene:",sum(unil==2)),
                     paste(">2 proteins/gene", sum(unil>2))), bty="n")
text(unil[unil>40], 100, labels=tgenes$name[unil>40], srt=90, pos=4, cex=.7)
dev.off()

## uniprot IDs per gene
unil <- unlist(lapply(strsplit(tgenes$uniprot, ";"),length))
plotdev(file.path(fig.path, "proteins_per_gene_uniprot"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unil, breaks=0:max(unil),
     xlab="uniprot IDs per `gene`",
     main=paste(sum(unil), "uniprot IDs"))
legend("topright", c(paste(" 1 protein /gene:",sum(unil==1)),
                     paste(" 2 proteins/gene:",sum(unil==2)),
                     paste(">2 proteins/gene", sum(unil>2))), bty="n")
text(unil[unil>40], 100, labels=tgenes$name[unil>40], srt=90, pos=4, cex=.7)
dev.off()

## isoforms per gene, e.g. multiple swissprot IDs
swissl <- unlist(lapply(strsplit(tgenes$swissprot, ";"),length))
plotdev(file.path(fig.path, "proteins_per_gene_swissprot"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(swissl, breaks=0:max(swissl),
     xlab="swissprot IDs per `gene`",
     main=paste(sum(swissl), "swissprot IDs"))
legend("topright", c(paste(" 1 protein /gene:",sum(swissl==1)),
                     paste(" 2 proteins/gene:",sum(swissl==2)),
                     paste(">2 proteins/gene", sum(swissl>2))), bty="n")
dev.off()
       
## isoforms per gene, e.g. multiple isoprot IDs
isol <- unlist(lapply(strsplit(tgenes$uniprot.isoforms, ";"),length))
plotdev(file.path(fig.path, "proteins_per_gene_uniprot_isoforms"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(isol, breaks=0:max(isol),
     xlab="uniprot isoform IDs per `gene`",
     main=paste(sum(isol), "uniprot isoform IDs"))
legend("topright", c(paste(" 1 protein /gene:",sum(isol==1)),
                     paste(" 2 proteins/gene:",sum(isol==2)),
                     paste(">2 proteins/gene", sum(isol>2))), bty="n")
#text(isol[isol>40], 100, labels=tgenes$name[isol>40], srt=90, pos=4, cex=.7)
dev.off()


## gene length
glen <- (abs(tgenes$end-tgenes$start)+1)
plotdev(file.path(fig.path, "gene_length"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(glen/1e3, xlab="`gene` length / kb", breaks=50,
     main=paste(nrow(tgenes),"protein_coding"))
text(glen[glen>1500000]/1e3, 1000,
     labels=tgenes$name[glen>1500000], srt=90, pos=4, cex=.7)
dev.off()

## go terms per gene
gol <- unlist(lapply(strsplit(tgenes$GO, ";"),length))
plotdev(file.path(fig.path, "go_per_gene"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(gol, xlab="GO terms per `gene`", breaks=seq(0,max(gol),1),
     main=paste(sum(gol),"GO terms"))
text(gol[gol>200], 100,
     labels=tgenes$name[gol>200], srt=90, pos=4, cex=.7)
dev.off()
## go slim terms per gene
gol <- unlist(lapply(strsplit(tgenes$GOslim, ";"),length))
plotdev(file.path(fig.path, "goslim_per_gene"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(gol, xlab="GOslim terms per `gene`", breaks=seq(0,max(gol),1),
     main=paste(sum(gol),"GOslim terms"))
text(gol[gol>40], 100,
     labels=tgenes$name[gol>40], srt=90, pos=4, cex=.7)
dev.off()
