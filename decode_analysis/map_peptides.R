
## MAP AAS TO ENSEMBL PROTEINS & TRANSCRIPTS
## and collect sequence & structure data at AAS.

## TODO: also load full list of unique SAAP/BP and do
## statistics on blast, AA and codon matching.

library(segmenTools)
options(stringsAsFactors=FALSE)


## DATA FROM  genomeBrowser, project folder data/mammary,
## run steps in data/mammary/setup.sh to create all data
## required here!

mam.path <- file.path(Sys.getenv("MAMDATA")) 
if ( mam.path=="" ) # author's local path
    mam.path <- "/home/raim/data/mammary"

feature.file <- file.path(mam.path,"features_GRCh38.110.tsv")

if ( !file.exists(feature.file) )
    stop("genome feature table file not found. This script requires ",
         "setup of genomic data via the genomeBrowser/data/mammary/setup.sh. ",
         "If you have setup this, please provide the path here as `mam.path` ",
         "and either change the path of the saap_mapped.tsv input or copy ",
         "it from processedData to additionalData.\n",
         "NOTE that the R analysis can still be run, since we provide ",
         "the output of this script, saap_mapped.tsv")


## DECODE DATA
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/decode"

dat.path <- file.path(proj.path,"originalData")
fig.path <- file.path(proj.path,"figures","saap_mapping")
out.path <- file.path(proj.path,"processedData")
out.file <- file.path(out.path,"saap_mapped.tsv")

## SET THIS VARIABLE TO GENERATE PDFs INSTEAD OF PNGs
ftyp <- "png"
##if ( !interactive() ) ftyp <- "pdf"

## SAAP/BP pairs
saapf <- file.path(out.path,"unique_saap.tsv")
## BP: protein blast results
bpmap <- file.path(out.path,"bp_mapped.tsv")

## protein fasta (Ensembl+Mutations, used for blast)
pfasta <- file.path(out.path,"all_proteins.fa")

## protein:transcript ID mapping
tpmap.file <- file.path(mam.path,"originalData","protein_transcript_map.tsv")

## coding region fasta
tfasta <- file.path(mam.path,"processedData","coding.fa")

## UTR Lengths
utr.file <- file.path(mam.path, "processedData", "transcript_utr_lengths.tsv")

## protein-transcript map
tpmap <- file.path(mam.path,"originalData","protein_transcript_map.tsv")
## CDS structure
cdsmap <- file.path(mam.path,"processedData","protein_cds_structure.dat")
cdspos <- file.path(mam.path,"processedData","protein_cds_coordinates.tsv")


## structure predictions
## s4pred
s4pred <- file.path(mam.path,"processedData",
                    "Homo_sapiens.GRCh38.pep.large_s4pred.fas.gz")
## iupred3/anchor2
iupred <- file.path(mam.path,"processedData","iupred3")

## uniprot<->ensembl mapping: required to find describePROT
uni2ens.file <- file.path(mam.path,"originalData","uniprot_ensembl.dat")

## describePROT
##  * sequence,
##  * MMseqs2 - Fast sequence alignment [PMID:30615063],
##  * ASAquick - Prediction of protein accessible surface area [PMID:27787824],
##  * DisoRDPbind - Prediction of disordered RNA, DNA, and protein
##    binding residues [PMID:26109352],
##  * SCRIBER - Prediction of protein binding residues [PMID:31510679],
##  * flDPnn - Prediction of intrinsically disordered residues [PMID:34290238].
descrp <- file.path(mam.path,"processedData","describePROT")

## pfam hits
pfam.file <- file.path(mam.path,"processedData",
                  "Homo_sapiens.GRCh38.pep.large_annotations.csv")
## pfam download from EBI
pfdl.file <- file.path(mam.path,"originalData", "9606.tsv.gz")
## pfam clans
clans.file <- file.path(mam.path,"originalData", "pfam", "Pfam-A.clans.tsv.gz")

## copied from the Biostrings package
GENETIC_CODE <- c(
    TTT="F", TTC="F",         
    TTA="L", TTG="L",         
    TCT="S", TCC="S", TCA="S", TCG="S", AGT="S", AGC="S",         
    TAT="Y", TAC="Y",         
    TAA="*", TAG="*", TGA="*",
    TGT="C", TGC="C",    
    TGG="W",         
    CTT="L", CTC="L", CTA="L", CTG="L",         
    CCT="P", CCC="P", CCA="P", CCG="P",         
    CAT="H", CAC="H",         
    CAA="Q", CAG="Q",         
    CGT="R", CGC="R", CGA="R", CGG="R", AGA="R", AGG="R",         
    ATT="I", ATC="I", ATA="I",
    ATG="M",         
    ACT="T", ACC="T", ACA="T", ACG="T",         
    AAT="N", AAC="N",         
    AAA="K", AAG="K",         
    GTT="V", GTC="V", GTA="V", GTG="V",         
    GCT="A", GCC="A", GCA="A", GCG="A",         
    GAT="D", GAC="D",         
    GAA="E", GAG="E",         
    GGT="G", GGC="G", GGA="G", GGG="G"
)

## figures
dir.create(fig.path, showWarnings=FALSE)



#### FIND POSITIONS of SAAP in protein and transcript

## GET ENSEMBL PROTEINS - from project mammary
fas <- readFASTA(pfasta, grepID=TRUE)

## get matching transcripts
trfas <- readFASTA(tfasta, grepID=TRUE)

## protein-transcript map
trmap <- read.delim(file=tpmap, header=FALSE, row.names=2)
## reverse map transcript-protein
pamrt <- matrix(rownames(trmap), ncol=1)
rownames(pamrt) <- trmap[,1]

## rename by protein names via trmap, for quick access
## of protein-specific transcript
names(trfas) <- pamrt[names(trfas),1]


## list of splice sites for each CDS
cds <- read.delim(cdsmap,header=FALSE, row.names=1)
cdl <- strsplit(cds[,1],";")
names(cdl) <- rownames(cds)
cdl <- lapply(cdl, as.numeric)
## genomic position of CDS
cpos <- read.delim(cdspos, header=FALSE, sep=" ")
colnames(cpos) <- c("ID","chr","start","end","strand")
posl <- split(cpos[,2:5], cpos[,1])
## NOTE: reverse order for minus strand
posl <- lapply(posl, function(x) x[order(x$start,
                                         decreasing=x$strand=="-"),])
## transcript positions: add UTRs
utr <- read.delim(utr.file, row.names=1)
rownames(utr) <- pamrt[rownames(utr),1]


## s4pred
s4p <- readFASTA(s4pred, grepID=TRUE)
## skip version number - TODO: use version numbers!
names(s4p) <- sub("\\.[0-9]+", "", names(s4p))
## reorder
cat(paste("s4pred:", sum(!names(fas)%in%names(s4p)), "proteins not found\n"))
s4p <- s4p[names(fas)]

## map iupred files here! 
iufiles <- list.files(pattern=paste0(".*iupred3.tsv.gz"), path=iupred)
names(iufiles) <- sub("\\..*","", iufiles)
iufiles <- iufiles[names(fas)]

cat(paste("iupred3", sum(is.na(iufiles)), "proteins not found\n"))

## read uniprot<->ensembl mapping for describPROT
uni2ens <- read.delim(uni2ens.file, header=FALSE)
uni2ens[,2] <- sub("\\..*", "", uni2ens[,2]) # remove ensembl version tag
uni2ens[,1] <- sub("-.*", "", uni2ens[,1]) # remove uniprot version tag
uni2ens <- uni2ens[uni2ens[,2]%in%names(fas),]

## PFAM 37.0 ANNOTATION from EBI
pfd <- read.delim(pfdl.file, skip=3, header=FALSE)
colnames(pfd) <- c("seq id", "alignment start",
                   "alignment end",
                   "FROM",#"envelope start",
                   "TO",#"envelope end",
                   "hmm acc", "hmm name", "type",
                   "hmm start", "hmm end", "hmm length",
                   "bit score", "E-value", "clan")
## add ensembl ID
pfd$ensembl <- uni2ens[match(pfd[,1], uni2ens[,1]),2]

## OWN PFAM/hmmer predictions - TODO, more likely to be correct!!
pfm <- read.csv(file=pfam.file, sep=";", fill=FALSE,
                header=FALSE, comment.char="#")
pfmh <- c(
    "target", "tid" , "tlen",            
    "query"   , "qid" , "qlen",                 
    "E-value"      , "score"     , "bias",                 
    "#"            , "of"        , "c-Evalue",             
    "i-Evalue"     , "score"     , "bias",            
    "from"         , "to"        , "from",                 
    "to"           , "FROM"      , "TO",                   
    "acc")#          , "description")
colnames(pfm) <- pfmh

## pfam clans: used to collapse domains below!
pclan <- read.delim(clans.file, row.names=1 , header=FALSE)
pfm$clan <- pclan[sub("\\..*","", pfm$tid),"V2"]
## same tag as EBI PFAM for No_clan
pfm$clan[pfm$clan==""] <- "No_clan"

## pfam w/o version number
pfm$pfam <- sub("\\..*","", pfm$tid)

## remove version number for access
pfm$ensembl <- sub("\\..*","", pfm$query)


### READ IN and MERGE unique SAAP/BP pairs and BP blast results
dat <- read.delim(saapf, header=FALSE)
colnames(dat) <- c("SAAP","BP")
bmap <- read.delim(bpmap)

dat <- merge(dat, bmap, by="BP", all=TRUE)


## result vectors
mut <- pos <- len <- cdn <- tps <- aaf <- aat <-  aas <-
    sss <- anc <- iup <- iubg <- anbg <-
        mmseq2 <- asaquick <- disordRDPbind <- scriber <- flDPnn <-
            pctx <- nctx <- pfam <- clan <- pfam.ebi <- clan.ebi <-
                rep(NA, nrow(dat))

## secondary structure frequencies in whole protein
sssbg <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(sssbg) <- c("C","E","H")

## genomic coordinate of AAS
gcoor <- matrix(NA, nrow=nrow(dat), ncol=3)
colnames(gcoor) <- c("chr","coor","strand")

## Quality Control
testit <- TRUE # TEST whether the mutation position is correct
use.regex <- TRUE #FALSE # use TRUE to test blast results vs. direct regex

## count errors
errtypes <- c("no protein", "using blast","no BP",
              "wrong s4pred len","wrong iupred3 len",
              "no describePROT",
              "wrong descPROT",
              "no transcript", "wrong codon",
              "AAS > CDS len")
wiup <- character()
errors <- matrix(0, nrow=nrow(dat), ncol=length(errtypes))
colnames(errors) <- errtypes
if ( !use.regex ) errors[,"using blast"] <- 1


## trouble shooting 20240424
pid="ENSP00000354876" # no transcript even though its in protein_transcript_map
## it's not in coding.fa, why??

## 20240807 #EGLELLK")
i=which(dat$BP=="ACQRPQLWQTIQTQGHFQLQLPPGK" &
        dat$SAAP=="ACQQPQLWQTIQTQGHFQLQLPPGK")

## 2002408

for ( i in 1:nrow(dat) ) {
   
    oid <- dat$protein[i] # ID with mutation index
    gid <- dat$ensembl[i] # original gene ID

    ## TODO: use MANE already here?
    ## problem is that we use oid for AAS mapping
    if ( FALSE ) {
        gid <- dat$MANE.protein[i] # original gene ID
        if ( is.na(gid)|gid=="" ) 
            gid <- dat$ensembl[i]
    }

    j <- which(names(fas)==oid)
    if ( length(j)==0 ) { # 
        cat(paste("WARNING:",i, oid, "no protein sequence found\n"))
        errors[i,"no protein"] <- 1
        next
    }
    if ( length(j)>1 ) { # atm not appearing
        stop(paste("WARNING:",i, oid, "more than one hit\n"))
        next
    }
    query <- unlist(dat[i,"BP"])
    target <- fas[[j]]$seq

    ## GET POSITION OF MUTATION
    ## 1: MAP PETPTIDE and test consistency with blast
    if ( use.regex ) {
        res <- gregexpr(query, target)
        if ( length(res)>1 ) {
            stop(paste("PROBLEM:", i, oid, 
                       "more than two hits in",j, names(fas)[j], "\n"))
        } else if ( res[[1]][1]==-1 ) {
            cat(paste("WARNING:",i, oid, "no match in", j, names(fas)[j],
                      ", using blast result with",
                      dat$mismatches[i]," mismatches, class:",dat$match[i],"\n"))
            errors[i,"using blast"] <- 1
            AAS <- dat[i,"sstart"]
        } else if ( res[[1]][1] != dat[i,"sstart"] ) {
            ## if this never occurs we could use blast results!
            stop(paste("PROBLEM:", i, oid, 
                       "regex doesn't match blast\n"))
        } else  AAS <- res[[1]][1]
    } else AAS <- dat[i,"sstart"]

    ## 2: AAS within protein
    saap <- unlist(dat[i,"SAAP"])
    mut[i] <- which(strsplit(saap,"")[[1]]!=strsplit(query,"")[[1]])
    pos[i] <- AAS + mut[i] -1  ## POSITION IN PROTEIN!
    len[i] <- length(unlist(strsplit(target,"")))

    ## GET AA SEQUENCE CONTEXT +/-25
    DST <- 25
    nsq <- nchar(target)
    rrng <- seq(-DST,DST,1)
    sq <- rep("-",length(rrng)) ## GAP
    names(sq) <- rrng
    ## GET RANGE AROUND AAS
    rng <- (pos[i]-DST):(pos[i]+DST)
    ## cut range to available
    rrng <- rrng[rng>0 & rng<=nsq]
    rng <- rng[rng>0 & rng<=nsq]
    sq[as.character(rrng)] <- unlist(strsplit(target,""))[rng]
    pctx[i] <- paste0(sq, collapse="")
   
    ## test mutation
    if ( testit ) { #
        if ( length(grep(query,target))==0 ) {
            cat(paste("WARNING:", i, gid, j, names(fas)[j],
                      "BP is not in target protein\n"))
            errors[i,"no BP"] <- 1
        } else {
            ntarg <- sub(query,saap,target)
            if ( strsplit(target,"")[[1]][pos[i]] ==
                 strsplit(ntarg,"")[[1]][pos[i]] ) {
                stop(paste("WARNING:", i, gid, j, names(fas)[j],
                           "no mutation detected\n"))
                
            }
        }
    }

    ## get PFAM CLANS - EBI download
    pfams <- pfd[which(pfd$ensembl==gid),,drop=FALSE]
    if ( nrow(pfams)> 0) {
        pidx <- which(pfams$FROM <= pos[i] & pfams$TO >= pos[i])
        if ( length(pidx)>0 ) {
            clan.ebi[i] <- paste0(unique(pfams[pidx, "clan"]), collapse=";")
            pfam.ebi[i] <- paste0(unique(pfams[pidx, "hmm acc"]), collapse=";")
        }
    }
    ## get PFAM CLANS - own hmmer prediction
    pfams <- pfm[which(pfm$ensembl==gid),,drop=FALSE]
    if ( nrow(pfams)> 0) {
        pidx <- which(pfams$FROM <= pos[i] & pfams$TO >= pos[i])
        if ( length(pidx)>0 ) {
            clan[i] <- paste0(unique(pfams[pidx, "clan"]), collapse=";")
            pfam[i] <- paste0(unique(pfams[pidx, "tid"]), collapse=";")
        }
    }

    ## GET s4pred
    s4s <- s4p[[gid]]
    if ( is.null(s4s) ) {
        cat(paste("WARNING:", i, "no s4pred found for", oid, "\n"))
        errors[,"wrong s4pred len"] <- "no s4pred result"
    } else {
        ## test length
        if ( nchar(s4s$seq)!=len[i] ) {
            cat(paste("WARNING:", i, "s4 prediction has wrong length for",
                      oid, "\n"))
            errors[,"wrong s4pred len"] <- 1
        } else {
            sss[i] <- substr(s4s$seq, pos[i], pos[i])
            ## background: complete protein sequence
            tb <- table(unlist(strsplit(s4s$seq,"")))
            sssbg[i,names(tb)] <- tb
        }
    }

    ## GET IUPRED3
    ## TODO: remove version tag from iupred files to avoid searching!
    ## or search once above
    iufile <- iufiles[gid]
    if ( is.na(iufile) )  {
        cat(paste("WARNING:", i, "no iupred3 found for", oid, "\n"))
        wiup <- c(wiup, paste(i, gid, "not found"))
        errors[i,"wrong iupred3 len"] <- "no iupred3 result"
    } else {
        iud <- read.delim(file.path(iupred,iufile), header=FALSE)
        if ( nrow(iud)!=len[i] ) {
            cat(paste("WARNING:", i, "iupred3 prediction has wrong length for",
                      oid, "iupred:",nrow(iud), "vs. protein:", len[i],"\n"))
            errors[i,"wrong iupred3 len"] <- 1
        } else {
            anc[i] <- iud[pos[i], 4] #anchor2
            iup[i] <- iud[pos[i], 3] #iupred3
            ## whole protein mean
            anbg[i] <- mean(iud[, 4])
            iubg[i] <- mean(iud[, 3])
        }

        ## NOTE: COPY TO SELECTED FOLDER FOR TRANSFER TO LAPTOP
        tof <- file.path(paste0(iupred,"_selected"),iufile)
        if ( Sys.info() ["nodename"]=="exon" ) 
            file.copy(from=file.path(iupred,iufile), to=tof, overwrite = FALSE)
    }

    ## get describePROT data
    if ( gid%in%uni2ens[,2] ) {
        idx <- which(uni2ens[,2]==gid)
        if ( length(idx)>1 ) {
            cat(paste("WARNING:", i, gid, "multiple uniprot hits,",
                      "taking first\n"))
            idx <- idx[1]
        }
        uid <- uni2ens[idx,1]
        ## open file
        dprt.file <- file.path(descrp, paste0(uid,".tsv.gz"))
        if ( !file.exists(dprt.file) ) {
                cat(paste("WARNING:", i,
                          "desribePROT missing for",
                          oid, uid, "\n"))
                errors[i,"no describePROT"] <- 1
        } else {
            dprt <- read.delim(dprt.file, header=FALSE)

            ## check sequence
            if ( nrow(dprt)!=len[i] ) {
                cat(paste("WARNING:", i,
                          "desribePROT has wrong length for",
                          oid, uid, "desribeP:",nrow(dprt), "vs. protein:",
                          len[i],"\n"))
                errors[i,"wrong descPROT"] <- 1
            }             
            mmseq2[i] <- dprt[pos[i], 2] 
            asaquick[i] <- dprt[pos[i], 3] 
            disordRDPbind[i] <- dprt[pos[i], 4] 
            scriber[i] <- dprt[pos[i], 5] 
            flDPnn[i] <- dprt[pos[i], 6] 
        }
    }
    
    ## GET CODON
    ## protein:transcript mapping
    ## load coding region get codon
    nt <- trfas[[gid]]
    if ( is.null(nt) ) {
        cat(paste("WARNING:", i, "no transcript found for", oid, "\n"))
        errors[i,"no transcript"] <- 1
        next
    }

    ## position of first codon position
    npos <- (pos[i]-1)*3+1
    ## retrieve codon from the transcript sequence
    codon <- substr(nt$seq, npos, npos+2)

    ## STORE FROM/TO AA
    aas[i] <- strsplit(target,"")[[1]][pos[i]]
    aaf[i] <- strsplit(query,"")[[1]][mut[i]]
    aat[i] <- strsplit(saap,"")[[1]][mut[i]]

    ## STORE CODON
    ## don't store wrong codon!
    if ( GENETIC_CODE[codon] != aaf[i] ) {
        cat(paste("WARNING:",i,"wrong codon", aaf[i],
                  "vs", codon, GENETIC_CODE[codon],"in", oid, "\n"))
        errors[i,"wrong codon"] <- 1
    } else {
        ## store codon
        cdn[i] <- codon
        ## store position in transcript
        if ( !gid%in%rownames(utr) )
            stop("missing UTR data")
        tps[i] <- npos + utr[gid,"start"] +1 # 2nd codon pos
    }
    
    ## GET NT SEQUENCE CONTEXT +/-25
    fsq <- nt$seq # transcript sequence
    nsq <- nchar(fsq)
    
    rrng <- seq(-DST*3,DST*3 +2,1)
    sq <- rep("-",length(rrng)) ## GAP
    names(sq) <- rrng

    ## GET RANGE AROUND AAS
    tpos <- pos[i]*3 -2
    rng <- (tpos-DST*3):(tpos+DST*3 +2)
    ## cut range to available
    rrng <- rrng[rng>0 & rng<=nsq]
    rng <- rng[rng>0 & rng<=nsq]
    sq[as.character(rrng)] <- unlist(strsplit(fsq,""))[rng]
    nctx[i] <- paste0(sq, collapse="")

    ## GET GENOMIC POSITION via CDS list and posl
    cds <- cdl[[gid]] # position of splice sites in transcript
    coors <- posl[[gid]] # genomic positions of CDS

    ## get second codon position
    npos <- npos+1

    exon <- which(cds>=npos)[1] # in which exon is the 2nd codon position
    if ( is.na(exon) ) {
        cat(paste("WARNING:",i,"wrong genomic position", npos, "vs", max(cds),
                  "in", gid, "\n"))
         errors[i,"AAS > CDS len"] <- 1
    } else {
    
        rpos <- npos # relative position in exon/CDS
        if ( exon>1 ) # subtract prior exons/CDS
            rpos <- npos-cds[exon-1]

        ## get genomic position of AAS 2nd codon position
        chr <- unique(coors$chr)
        strand <- unique(coors$strand)
        
        if ( strand=="-") {
            coor <- coors[exon, "end"] - rpos +1
        } else {
            coor <- coors[exon, "start"] + rpos -1
        }
        gcoor[i,] <- c(chr, coor, strand)
    }
     
    ## report
    if ( !interactive() ) {        
        tcodons <- paste(names(which(GENETIC_CODE%in%aat[i])),collapse=";")
        cat(paste("DONE:", i, codon, aaf[i], GENETIC_CODE[codon], "->", aat[i],
                  tcodons,  "\n"))
    }
    
    if ( aas[i] != aaf[i] ) 
        warning(paste("WARNING: BP pos from AA:", aaf[i], "vs",
                  aas[i], GENETIC_CODE[codon], gid, "\n"))
}






## relative position
rpos <- pos/len
## bind to data frame
bgsss <- sssbg
colnames(bgsss) <- paste0(colnames(bgsss), ".protein")
dat <- cbind(dat,
             site=mut,  # mutated site within peptide
             pos=pos,   # position of AAS within protein
             len=len,
             rpos=pos/len,
             from=aaf, to=aat,
             codon=cdn, 
             tpos=tps,   # position of AAS codon within transcript
             gcoor,      # genome coordinate of AAS
             s4pred=sss, bgsss, iupred3=iup, iupred3.protein=iubg,
             anchor2=anc, anchor2.protein=anbg,
             MMSeq2=mmseq2,
             ASAquick=asaquick,
             DisoRDPbind=disordRDPbind,
             SCRIBER=scriber,
             flDPnn=flDPnn,
             AA=pctx,
             NT=nctx,
             pfam=pfam,
             clan=clan,
             pfam.ebi=pfam.ebi,
             clan.ebi=clan.ebi)

## FILTER DATA

## missing position - no match found from peptide in protein!!
## TODO: perhaps test LONGEST of the leading razor proteins.
rm <- which(is.na(dat$pos)) #dat$pos<0)
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    ##errors <- errors[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP without protein match\n"))
}

## remove where length or position are NA - same as "no match"!
rm <- which(is.na(dat$len))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    ##errors <- errors[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP no protein length\n"))
}

## remove where no codon was assigned
## no transcript found or wrong codon
rm <- which(is.na(dat$codon))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP w/o codon\n"))
}

## remove where no codon was assigned
## no transcript found or wrong codon
rm <- which(is.na(dat$chr))
if ( length(rm)>0 ) {
    ##dat <- dat[-rm,]
    cat(paste("WARNING:", length(rm), "SAAP w/o genome coordinates\n"))
}

### ERROR AND RECOVERY STATS

## NOTE:
## * "no BP" and "using blast" are the same, and all are "bad" or "wrong"
##   blast matches,
## * TCHVEHPSLQNPITVEW:TCHVEHPSLQSPITVEW:ENSP00000382018
##   is classified as "wrong" in blast, because it is associated with
##   a scaffold gene (is.na(all$gene)) and not with the official genome gff3.

errs <- errors#[,apply(errors,2,sum)>0]
er <- c("unique BP/SAAP"=nrow(dat),
        "unique SAAP"=sum(!duplicated(dat$SAAP)),
        "unique BP"=sum(!duplicated(dat$BP)),
        table(dat$match)[c("good","bad","wrong")],
        "no gene"=sum(dat$nogene, na.rm=TRUE),
        "with mutation"=sum(dat[,"protein"]!=dat[,"ensembl"],na.rm=TRUE),
        apply(errs,2,sum))
plotdev(file.path(fig.path,"mapping_errors"),
        res=300, width=4, height=3.5, type=ftyp)
par(mai=c(1.5,.5,.5,.1), mgp=c(1.3,.3,0), tcl=-.25)
bp <- barplot(er,las=2)
text(bp+.25, er+.05*diff(par("usr")[3:4]),
     labels=er, pos=3, xpd=TRUE, cex=.8, srt=90)
dev.off()

## TODO: better characterize matches to mutated proteins
recovery <- c("SAAP"=nrow(dat),
              "protein"=sum(!is.na(dat[,"protein"])),
              ##"length"=sum(!is.na(dat[,"len"])),
              "s4pred"=sum(!is.na(dat[,"s4pred"])),
              "iupred3"=sum(!is.na(dat[,"iupred3"])),
              "genome"=sum(!is.na(dat[,"chr"])),
              "codon"=sum(!is.na(dat[,"codon"])))
plotdev(file.path(fig.path,"mapping_recovery"),
        res=300, width=4, height=3.5, type=ftyp)
par(mai=c(.75,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25)
bp <- barplot(recovery,las=2,ylab="")
text(bp, recovery, recovery, pos=3, xpd=TRUE, cex=.8)
dev.off()

types <- c("SAAP"=nrow(dat),
           "MANE"=sum(!is.na(dat[,"MANE.protein"])),
           "extracellular"=sum(dat[,"extracellular"], na.rm=TRUE),
           "IG"=sum(dat[,"IG"], na.rm=TRUE),
           "albumin"=sum(dat[,"albumin"], na.rm=TRUE),
           "globin"=sum(dat[,"globin"], na.rm=TRUE))

types <- sort(types, decreasing=TRUE)
plotdev(file.path(fig.path,"mapping_types"),
        res=300, width=4, height=3.5, type=ftyp)
par(mai=c(.75,.5,.25,.1), mgp=c(1.3,.3,0), tcl=-.25, cex=.8)
bp <- barplot(types,las=2,ylab="")
text(bp, types, types, pos=3, xpd=TRUE)
dev.off()

## inspect some overlaps
table(dat[,"extracellular"], dat[,"globin"])
table(dat[,"extracellular"], dat[,"albumin"])
table(dat[,"extracellular"], dat[,"IG"])


### WRITE OUT TABLE with positions for downstream analysis
if ( !interactive() ) {
    write.table(dat, file=out.file,
                sep="\t", quote=FALSE, na="", row.names=FALSE)
}

