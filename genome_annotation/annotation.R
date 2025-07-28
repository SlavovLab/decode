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
## T<-P map
tpmap.file <- file.path(mam.path,"originalData", "protein_transcript_map.tsv")

## tRNAs - TODO:  integrate with codon:anticodon mapping
trna.file <- file.path(mam.path,"originalData",
                       "hg38-tRNAs.bed")

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
## other data
ortho.file <- file.path(mam.path,"processedData","human_h38_20231204.tsv")
tfeb.file <- file.path(mam.path, "originalData", "sardiello09_tabS3.txt")

### TODO: EXPERIMENTAL DATA

## human protein halflives - @Mathieson2018
hhlf.file <- file.path(mam.path,"originalData",
                       "41467_2018_3106_MOESM5_ESM.xlsx")
## protein lengths
hlen.file <- file.path(mam.path,"processedData",
                       "protein_length.tsv")

## @Yang2022: thermal stability prediction
## https://structure-next.med.lu.se/ProTstab2/
protstab.file <- file.path(mam.path,"originalData",
                           "ProTstab2_human.csv")

## @Savitski2014: Tracking cancer drugs in living cells by thermal
## profiling of the proteome
thermo.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS11.xlsx")
## Diff melting point with ATP
thatp.file <- file.path(mam.path, "originalData",
                         "savitski14_tableS3.xlsx")

## 20S targets - @Pepelnjak2024
pepe24.file <- file.path(mam.path,"originalData",
                         "44320_2024_15_moesm1_esm.xlsx")



## output: mapping statistics
feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

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


### ORTHOLOGS
ortho <- read.delim(ortho.file)

## no duplicate yeast genes
sum(duplicated(ortho$Orthologs))

## split orthologs into vectors
orthl <- lapply(strsplit(ortho$human_h38_20231204,","), trimws)
names(orthl) <- ortho$Orthologs
## expand list to full map
orthm <- list2df(orthl)
## replace ensembl version number!!
orthm[,2] <- sub("\\.[0-9]+","",orthm[,2])


## since there are no duplicated human orthologs
## we can convert it to a named vector
sum(duplicated(orthm[,2]))
orthv <- orthm[,1]
names(orthv) <- orthm[,2]

## map to protein list!
orthl <- lapply(prol, function(x) unique(orthv[x]))
orthl <- lapply(lapply(orthl, function(x) unlist(strsplit(x, ","))), trimws)
## rm NA
orthl <- lapply(orthl, function(x) x[!is.na(x)])


## check a gene with several yeast orthologs
if ( interactive() ) {
    id="ENSG00000292334"
    orthl[id]
}

cat(paste("found", length(unlist(orthl)), "yeast orthologs\n"))

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

#### COLLECT EXPERIMENTAL DATA
## 



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
yeast <- unlist(lapply(orthl[match(genes$ID, names(orthl))], collaps))

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
                GOslim=go, GO=gof,
                yeast=yeast)


## add tRNA annotation from gtRNADB
## NOTE: bed - 0-based starts
trna <- read.delim(trna.file, header=FALSE)
colnames(trna)[1:6] <- c("chr","start","end", "ID", "V5","strand")
trna$start <- trna$start+1 # correct 0-based coor start

## inspect anticodons
table(sub("tRNA-Thr-","",sub("-[0-9]*-[0-9]*$","",grep("Thr", trna$ID, value=TRUE))))
## NOTE: no ACC anticodon (GGT), even though it is the most frequent codon,
## likely: tRNA-Thr-AGT is converted to IGT, and I can pair with U,C,A

## NOTE: tRNA-Val-AAC could decode thr
## see Table 1. Suggested Codon Usage Preferences in Each Tissue Relative to Brain According to tRNA Expression of https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020221

## add manual annotations


                                              
write.table(agenes,
            file=feature.file,
            sep="\t", row.names=FALSE, quote=FALSE, na="")


### STATISTICS and PLOTS
## statistics table and figures
## transcripts, proteins, UTRs per gene etc.
## TODO: yeast orthologs per human gene and reverse


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
hist(unlist(tesl)/1e3, xlab="TES distance from `gene` end / kb", xlim=c(-100,0),
     breaks=1000, main="transcript end sites")
dev.off()
plotdev(file.path(fig.path, "TSS_distance_zoom"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(unlist(tssl)/1e3, xlab="TSS distance from `gene` start / kb", xlim=c(0,100),
     breaks=1000, main="transcript start sites")
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

## TODO:
## * gene length vs. all its transcript lengths,
## * protein length vs transcript lengths (exons).

## generate markdown table of full feature table types
ttab <- sort(table(genes$type), decreasing=TRUE)
## add dummy line
ttab <- c(ttab, NA)
cat(paste(paste("|",names(ttab)[1:20], "|", ttab[1:20],"||",
                names(ttab)[21:40], "|", ttab[21:40],"|"),collapse="\n"))

## example for a 1Mb gene with long distance TSS/TES - many associated proteins
## Gene: CAMTA1 - calmodulin binding transcription activator 1
## http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000171735;r=1:6785454-7769706;t=ENST00000700415

## most isoforms
## Gene: ANK2 - 

## example of a gene with many isoforms
## Gene: SGIP1 - SH3GL interacting endocytic adaptor 1
## https://www.uniprot.org/uniprotkb/A0A804HHZ7/entry

### FUNCTIONAL ANALYSIS

## get GOslim table for use with clusterAnnotation
got <- parseAnnotationList(tgenes[,c("ID","GOslim")]) 
## replace GO IDs by terms
trms <- terms[,2]
names(trms) <- terms[,1]
colnames(got) <- trms[colnames(got)]


## load yeast clustering and do GO analysis for orthologs
ygene.file <- "~/data/yeast/feature_R64-1-1_20110208_withclusters.csv"
ygenes <- read.delim(ygene.file)
ygenes <- ygenes[ygenes$type=="gene",]

## get yeast gene<>cluster map
ycls <- ygenes$CL_rdx
names(ycls) <- ygenes$ID


## TODO: handle case of multiple yeast genes!!
## -> clusterAnnotation for multi-multi mapping instead!!
## simpler: use list, get unique, and classify multi as multi,
## or add as extra class
## FOR NOW: just take the first

## map to ortholog column!
ygns <- strsplit(tgenes$yeast, ";")
ygns <- lapply(ygns, function(x) ycls[x])
names(ygns) <- tgenes$ID

hcls <- unlist(lapply(ygns, function(x) x[1]))
names(hcls) <- tgenes$ID
##hcls <- hcls[rownames(got)]
hcls[is.na(hcls)] <- "na"

## cluster sorting
cls.srt <- sort(unique(hcls))
cls.srt.major <- cls.srt[grep("[A-Z]", cls.srt)]
cls.srt.minor <- cls.srt[!cls.srt%in%cls.srt.major]
cls.nms <- cls.srt
names(cls.nms) <- cls.srt
name.map <- c(A="Ribi",
              AB="RP",
              B="AA",
              B.C="TCA+S",
              B.D="AA/N",
              C="mRP",
              D="S/C")
cls.nms[names(name.map)] <- name.map

## calculate enrichment!
ovl <- clusterAnnotation(cls=hcls, data=got, cls.srt=cls.srt)

## cut at minimal level of significance
ovc <- sortOverlaps(ovl, p.min=1e-5, cut=TRUE)

plotdev(file.path(fig.path,"yeast_orthologs_goslim"),
        typ=ftyp, res=res, width=10, height=15)
par(mai=c(.7,4,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovc$p.value), labels=cls.nms[colnames(ovc$p.value)],
    las=2)
dev.off()

## major clusters only
ovm <- sortOverlaps(ovl, axis=1, srt=cls.srt.major)
ovm <- sortOverlaps(ovm, p.min=1e-5, cut=TRUE)

plotdev(file.path(fig.path,"yeast_orthologs_goslim_major"),
        typ=ftyp, res=res, width=7, height=8)
par(mai=c(.7,4,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovm, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovm$p.value), labels=cls.nms[colnames(ovm$p.value)],
    las=2)
dev.off()

## minor clusters only
ovm <- sortOverlaps(ovl, axis=1, srt=cls.srt.minor)
ovm <- sortOverlaps(ovm, p.min=1e-5, cut=TRUE)

plotdev(file.path(fig.path,"yeast_orthologs_goslim_minor"),
        typ=ftyp, res=res, width=7, height=12)
par(mai=c(.7,4,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovm, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovm$p.value), labels=cls.nms[colnames(ovm$p.value)],
    las=2)
dev.off()

## load TFEB target genes

tfeb <- read.delim(tfeb.file, header=FALSE)[,1]

## TODO: some are missing!
tfeb[!tfeb%in%tgenes$name]


tidx <- c(grep("lysoso", terms[,2]),
          grep("autophagy", terms[,2]),
          grep("programmed cell death", terms[,2]))

## matrix: experimental target genes, and GO
lysm <- matrix(FALSE, ncol=1+length(tidx), nrow=nrow(tgenes))
lysm[,1] <- tgenes$name%in%tfeb
for ( j in seq_along(tidx) )
    lysm[grep(tidx[j], tgenes$GOslim),j+1] <- TRUE
colnames(lysm) <- c("TFEB",terms[tidx, 2])

## calculate enrichment!
ovl.tfeb <- clusterAnnotation(cls=hcls, data=lysm, cls.srt=cls.srt)

## cut and sort along clustering
ovs <- sortOverlaps(ovl.tfeb, axis=1, srt=cls.srt.major)
#ovs <- sortOverlaps(ovs, p.min=1e-5, cut=TRUE)
#ovs <- sortOverlaps(ovs, p.min=1e-7)

plotdev(file.path(fig.path,"yeast_orthologs_tfeb"),
        typ=ftyp, res=res, width=7, height=1.4+.25*ncol(lysm))
par(mai=c(.7,2,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl.tfeb, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovl.tfeb$p.value),
     labels=cls.nms[colnames(ovl.tfeb$p.value)],
    las=2)
dev.off()

plotdev(file.path(fig.path,"yeast_orthologs_tfeb_major"),
        typ=ftyp, res=res, width=7, height=1.4+.25*ncol(lysm))
par(mai=c(.7,4,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovs$p.value), labels=cls.nms[colnames(ovs$p.value)],
    las=2)
dev.off()

## MOST TFEB target genes are not in ortholog list
table(tgenes$yeast!="", lysm[,1])


### LOAD clustering of MCF10A time series

gm.file <- file.path(mam.path, "processedData", "gutierrezmonreal16",
                     "gutierrezmonreal16_clustering.rda")

## TODO: nicer clustering files for gm16 data
tmp <- cls.srt
load(gm.file)
cls.srt <- tmp
cset1 <- cset

## manual sorting based on overlap
cls1.srt <- as.character(c(11,4,2,9,13,
                           8,5,1,3,7,12,10,6))
cls7.srt <- as.character(c(1,5,9,10,8,
                           2,7,11,3,
                           6,4))

## map MC10A clustering to tgenes
for ( ct in c(1,7) ) {

    cset <- get(paste0("cset",ct))
    gm16.cls.srt <-  get(paste0("cls",ct,".srt"))
    
    gm16.cls <- cset$clusters[,selected(cset)]
    gm16.cls <- gm16.cls[match(tgenes$name, names(gm16.cls))]
    gm16.cls[is.na(gm16.cls)] <- "na"


    ## calculate enrichment!
    ovl.gm16 <- clusterCluster(target=hcls, query=gm16.cls, t.srt=cls.srt,
                               q.srt=gm16.cls.srt)
    
    ## cut and sort along clustering
    ovs <- sortOverlaps(ovl.gm16, axis=1, srt=cls.srt.major)
    
    
    plotdev(file.path(fig.path,paste0("yeast_orthologs_gm16_",ct)),
            typ=ftyp, res=res, width=7, height=5)
    par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovl.gm16, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                 xlab=NA, ylab=NA, axis=2)
    imap(axis, 1, at=1:ncol(ovl.gm16$p.value),
         labels=cls.nms[colnames(ovl.gm16$p.value)],
         las=2)
    dev.off()

    plotdev(file.path(fig.path,paste0("yeast_orthologs_gm16_",ct,"_major")),
            typ=ftyp, res=res, width=4, height=5)
    par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
    plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
                 xlab=NA, ylab=NA, axis=2)
    imap(axis, 1, at=1:ncol(ovs$p.value), labels=cls.nms[colnames(ovs$p.value)],
         las=2)
    dev.off()
}


## * cell cycle genes, bin by phase

cdc.file <- file.path(mam.path, "originalData", "human_periodic.tsv")
cdc <- read.delim(cdc.file)


## map to tgenes
prol <- strsplit(tgenes$proteins, ";")
names(prol) <- tgenes$ID
g2p <- list2df(prol)
rownames(g2p) <- g2p$vectors
cdc.genes <- g2p[cdc$gene,1]
cdc <- cdc[match(tgenes$ID, cdc.genes),]

plotdev(file.path(fig.path, "cdc_peaktime"),
        width=W, height=H, type=ftyp, res=res)
par(mai=c(.5,.5,.15,.1), tcl=-.25, mgp=c(1.3,.3,0))
hist(cdc$peaktime, xlab="cell cycle peak time",
     main=paste(sum(!is.na(cdc$peaktime)), "CDC genes at `cyclebase`"))
dev.off()

cdc.cls <- as.character(cut(cdc$peaktime+1, breaks=seq(0,100,20)))
cdc.cls[is.na(cdc.cls)] <- "na"
cdc.cls.srt <- c(levels(cut(cdc$peaktime+1, breaks=seq(0,100,20))),"na")

## calculate enrichment!
ovl.cdc <- clusterCluster(target=hcls, query=cdc.cls, t.srt=cls.srt,
                          q.srt=cdc.cls.srt)
## cut and sort along clustering
ovs <- sortOverlaps(ovl.cdc, axis=1, srt=cls.srt.major)


plotdev(file.path(fig.path,paste0("yeast_orthologs_cdc")),
        typ=ftyp, res=res, width=7, height=5)
par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovl.cdc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovl.cdc$p.value),
     labels=cls.nms[colnames(ovl.cdc$p.value)],
     las=2)
dev.off()

plotdev(file.path(fig.path,paste0("yeast_orthologs_cdc_major")),
        typ=ftyp, res=res, width=4, height=5)
par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovs$p.value), labels=cls.nms[colnames(ovs$p.value)],
     las=2)
dev.off()

#### YEAST GENES - OSCI vs. CDC
## ignoring human: TODO: move to yeast project,
## and align with clusterTranscriptomes.R

## TODO: check cyclebase
ycdc.file <- file.path(mam.path, "originalData", "cerevisiae_periodic.tsv")
ycdc <- read.delim(ycdc.file)

## map to yeast genes
ycdc <- ycdc[match(names(ycls), ycdc$gene),]

cdc.cls <- as.character(cut(ycdc$peaktime+1, breaks=seq(0,100,20)))
cdc.cls[is.na(cdc.cls)] <- "na"
cdc.cls.srt <- c(levels(cut(ycdc$peaktime+1, breaks=seq(0,100,20))),"na")

## calculate enrichment!
ovl.cdc <- clusterCluster(target=ycls, query=cdc.cls, t.srt=cls.srt,
                          q.srt=cdc.cls.srt)
## cut and sort along clustering
ovs <- sortOverlaps(ovl.cdc, axis=1, srt=cls.srt.major)


plotdev(file.path(fig.path,paste0("yeast_cdc")),
        typ=ftyp, res=res, width=7, height=5)
##par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
par(mai=c(.9,.9,.6,.6), mgp=c(1.3,.4,0),tcl=-.25)
plotOverlaps(ovl.cdc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovl.cdc$p.value),
     labels=cls.nms[colnames(ovl.cdc$p.value)],
     las=2)
mtext("cyclebase, phase bins", 2, 3.3)
dev.off()

plotdev(file.path(fig.path,paste0("yeast_cdc_major")),
        typ=ftyp, res=res, width=3.5, height=3.5)
par(mai=c(.7,.9,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovs$p.value), labels=cls.nms[colnames(ovs$p.value)],
     las=2)
mtext("cyclebase, phase bins", 2, 3.3)
dev.off()


## load and map new clusters
ynew.file <- "~/work/ChemostatData/figures/pwmavg/genedata.tsv"
ynew <- read.delim(ynew.file)

new.cls <- ynew$CL.segment[match(ygenes$ID, ynew$ID)]
new.cls[is.na(new.cls)] <- "na"

new.srt <- c(as.character(c(1:10,0)),"na")
new.nms <- new.srt
names(new.nms) <- new.srt
new.nms["1"] <- "Ribi"
new.nms["2"] <- "RP"
new.nms["3"] <- "AA"
new.nms["7"] <- "mRP"
new.nms["8"] <- "S/C"
new.nms["4"] <- "nc1"
new.nms["5"] <- "nc2"
new.nms["6"] <- "nc3"
new.nms["9"] <- "nc4"
new.nms["10"] <- "nc5"
new.nms["0"] <- "nc6"

new.srt.major <-  new.srt[grep("^n", new.nms[new.srt], invert=TRUE)]

table(new.nms[new.cls])

## calculate enrichment!
ovl.cdc <- clusterCluster(target=new.cls, query=cdc.cls, t.srt=new.srt,
                          q.srt=cdc.cls.srt)
## cut and sort along clustering
ovs <- sortOverlaps(ovl.cdc, axis=1, srt=new.srt.major)


plotdev(file.path(fig.path,paste0("yeast_cdc_new")),
        typ=ftyp, res=res, width=7, height=5)
##par(mai=c(.7,.7,.7,.7), mgp=c(1.3,.3,0), tcl=-.25)
par(mai=c(.9,.9,.6,.6), mgp=c(1.3,.4,0),tcl=-.25)
plotOverlaps(ovl.cdc, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovl.cdc$p.value),
     labels=new.nms[colnames(ovl.cdc$p.value)],
     las=2)
mtext("cyclebase, phase bins", 2, 3.3)
dev.off()

plotdev(file.path(fig.path,paste0("yeast_cdc_major_new")),
        typ=ftyp, res=res, width=3.5, height=3.5)
par(mai=c(.7,.9,.5,.5), mgp=c(1.3,.3,0), tcl=-.25)
plotOverlaps(ovs, p.min=1e-10, p.txt=1e-5, show.total=TRUE,
             xlab=NA, ylab=NA, axis=2)
imap(axis, 1, at=1:ncol(ovs$p.value), labels=new.nms[colnames(ovs$p.value)],
     las=2)
mtext("cyclebase, phase bins", 2, 3.3)
dev.off()


## TODO: compare cyclebase genes via orthologs
