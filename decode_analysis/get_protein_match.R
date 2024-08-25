
## load blast of base peptides against ensembl+mutationes, and
## find best matching protein for each unique BP

library(segmenTools)
options(stringsAsFactors=FALSE)



## DECODE DATA
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/decode"

dat.path <- file.path(proj.path,"originalData")
bp.file <- file.path(proj.path,"processedData","unique_bp_blast.tsv")
out.file <- file.path(proj.path,"processedData","bp_mapped.tsv")

## DATA FROM  genomeBrowser, project folder data/mammary,
## run steps in data/mammary/setup.sh to create all data
## required here!
##mam.path <- "/home/raim/data/mammary"
##feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")
##tpmap.file <- file.path(mam.path,"originalData","protein_transcript_map.tsv")

## USING LOCAL COPIES INSTEAD (for publication git)!
feature.file <- file.path(proj.path,"additionalData","features_GRCh38.110.tsv.gz")
tpmap.file <- file.path(proj.path,"additionalData","protein_transcript_map.tsv.gz")


### LOAD & ANNOTATE BLAST RESULTS


if ( !file.exists(bp.file) )
    stop("input file is not available; make sure that either DECODE is set in",
         "the calling environment or that the variable proj.path in this",
         "script")


bp <- read.delim(bp.file, header=FALSE)
colnames(bp) <- c("BP","protein","identity", "mismatches",
                  "alen", "qlen", "slen", "sstart", "send", "e", "bitscore")
bp$ensembl <- sub("_.*", "", bp$protein)

## load feature file
features <- read.delim(feature.file)
## use number of GO annotations as a filter criterium: best annotated
gol <- unlist(lapply(strsplit(features$GO, ";"),length))

## load protein-transcript map
tpmap <- read.delim(file=tpmap.file, header=FALSE, row.names=2)
## reverse map transcript-protein
pampt <- matrix(rownames(tpmap), ncol=1)
rownames(pampt) <- tpmap[,1]

## map proteins to feature table!
plst <- strsplit(features$proteins, ";")
gidx <- rep(1:nrow(features), lengths(plst))
names(gidx) <- unlist(plst)

## genes and transcripts for each protein
bp$gene <- features$ID[gidx[bp$ensembl]]
bp$transcript <- tpmap[bp$ensembl,1]

## ADD MANE transcript/protein
## MANE: Matched Annotation from NCBI and EBI
## see ensembl MANE project: consistently annotated
tmp <- features$MANE[gidx[bp$ensembl]]
tmp[tmp==""] <- NA
bp$MANE <- features$MANE[gidx[bp$ensembl]] ##unlist(tmp)
bp$MANE.protein <- pampt[match(bp$MANE,rownames(pampt)),1]
bp$numGO <- gol[gidx[bp$ensembl]]


### SORT & FILTER BLAST HITS

## sort by quality - such that the best option is always
## selected (e.g. protein==MANE) when filtering duplicates
## in this order:
## identity > mismatch>
## MANE present > initial match==MANE >
## shortest> number of annotations
bp <- bp[order(bp$identity,
               -bp$mismatches,
               !is.na(bp$MANE.protein),
               bp$protein==bp$MANE.protein,
               bp$slen,
               bp$numGO, decreasing=TRUE),]

## FIRST HANDLE PERFECT HITS 

all <- bp ## store to analyze missing below

## filter - only those where we found a gene, others
## w/o gene are encoded on scaffolds, TODO: check all.
ina <- is.na(bp$gene)
cat(paste("removing", sum(ina), "blast hits without a mapped gene!\n"))
bp <- bp[!ina,] # proteins annotated to a scaffold


mm <- bp$mismatches>0 | bp$identity <100
cat(paste("removing", sum(mm), "hits with mismatches or <100% identity\n"))
bp <- bp[!mm,]

## THE WINNER TAKES IT ALL/RICH GET RICHER APPROACH.
## TODO: check this approach, perhaps there is a better way?

## number of AAS per MANE protein
## this is used to select the "winner" protein match below.
app <- table(bp$MANE.protein) 

## split into lists by BP
bpl <- split(bp, paste(bp$BP)) 

## for each BP select the MANE protein with the
## the most global matches.
winner <- list()
for ( i in seq_along(bpl) ) {
    x <- bpl[[i]]
    mane <- x$MANE.protein
    if ( length(unique(mane))>1 ) {
        ## only keep the mane protein with most hits globally
        cat(paste(i, x$BP[1],
                  "selecting MANE protein from multiple hits",
                  nrow(x)))
        manemax <- which(mane==mane[which.max(app[mane])])
        x <- x[manemax,,drop=FALSE] 
        cat(paste("->", nrow(x), "\n"))
    }
    winner[[i]] <- x[1,,drop=FALSE] ## ONLY TAKE FIRST
}
best <- do.call(rbind, winner)

## HANDLE REST
rest <- all[!all$BP%in%c(best$BP),]

## SIMPLY REMOVE DUPLICATE BP and rely on above sorting
dupb <- duplicated(paste0(rest$BP))
cat(paste("removing", sum(dupb), "multiple matches from rest\n"))
rest <- rest[!dupb,]

## combine again and tag quality
all <- rbind(cbind(best, match="good"),
             cbind(rest, match="bad"))

## tag really bad hits
all$match[all$mismatch>=3] <- "wrong"

## tag no gene
all$nogene <- FALSE
all$nogene[is.na(all$gene)] <- TRUE



### TAG CATEGORIES

type <- features$type[gidx[all$ensembl]]
all$IG <- FALSE
all$IG[grep("^IG",type)] <- TRUE

desc <- features$description[gidx[all$ensembl]]
all$albumin <- FALSE
all$albumin[grep("albumin",desc)] <- TRUE
all$globin <- FALSE
all$globin[grep("globin",desc)] <- TRUE

## load GO terms
## extracellular region: GO:0005576
## |- extracellular space: GO:0005615
gots <- features$GO[gidx[all$ensembl]]
all$extracellular <- FALSE
all$extracellular[grep("GO:0005576", gots)] <- TRUE

## exclude tag  - TODO: any K/R                 
all$exclude <- all$globin|all$albumin|all$IG

## ADD GENE NAME
all$name <- features$name[gidx[all$ensembl]]

if ( !interactive() ) {
    write.table(all, file=out.file, sep="\t",
                quote=FALSE, na="", row.names=FALSE)
}

