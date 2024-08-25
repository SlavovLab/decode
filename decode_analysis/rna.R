
## mRNA MODIFICATIONS at AMINO ACID SUBSTITUTION SITES  

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")

## local output path
rfig.path <- file.path(fig.path,"rna")
dir.create(rfig.path, showWarnings=FALSE)

## use only unique AAS sites where we have a codon
usite <- asite[!is.na(asite$codon),]

## use only AAS mapped to a MANE transcripts
usite <- usite[usite$transcript%in%genes$MANE,]

## prepare AAS data and coordinates


## get chromosome index
chrMap <- read.delim(chr.file)
chrS <- c(0,cumsum(as.numeric(chrMap[,3])))
chrIdx <- chrMap[,1]
names(chrIdx) <- chrMap[,2]

## convert AAS chromosome to index
usite$chr <- chrIdx[usite$chr]

## AAS 2nd codon position
bd2dx <- coor2index(usite[,c("chr","coor")], chrS=chrS)[,2]
## AAS 1st codon position
bd1 <- cbind(chr=usite[,c("chr")],
             coor=usite[,c("coor")] + ifelse(usite[,"strand"]=="+",-1,1))
bd1dx <- coor2index(bd1, chrS=chrS)[,2]
## AAS 3rd codon position
bd3 <- cbind(chr=usite[,c("chr")],
             coor=usite[,c("coor")] - ifelse(usite[,"strand"]=="+",-1,1))
bd3dx <- coor2index(bd3, chrS=chrS)[,2]

### COLLECT VARIOUS PSI DATA

PSI <- NULL

### RNA pseudouridylation data sent by
### Oleksandra Fanari <fanari.o@northeastern.edu>
psi <- NULL
for ( i in excel_sheets(psi.file) ) {
    
    tmp <-as.data.frame(read_xlsx(psi.file, sheet=i))

    ## collapse duplicate sites
    ## check whether
    upos <- paste(tmp$chr, tmp$position)
    dpos <-upos[which(duplicated(upos))]

    ## all annotations as list
    danl <- sapply(dpos, function(x) tmp$Annotation[which(upos==x)])
    ## filter those that are present in protein name vector
    danl <- lapply(danl, function(x) x[x%in%pnms])

    ## TODO: further reduce to a single annotation
    ## e.g. by length

    dann <- sapply(dpos, function(x)
        paste(tmp$Annotation[which(upos==x)], collapse=";"))
    fpos <- sapply(dpos, function(x)
        which(upos==x)[1])

    ## add annotation lists to be filtered later
    tmp$Annotation.all <- tmp$Annotation.best <- tmp$Annotation.first <-
        tmp$Annotation
    tmp$Annotation.all[fpos] <- dann
    tmp$Annotation.best[fpos] <- unlist(lapply(danl, paste, collapse=";"))
    tmp$Annotation.first[fpos] <- unlist(lapply(danl, function(x) x[1]))

    ## remove all duplicated
    tmp <- tmp[!duplicated(upos),]
    
    psi <- rbind(psi, cbind(tmp, cell=i))
}
colnames(psi) <- sub("position","coor", colnames(psi))

## convert chromosomes ID to index via segmenTools'
## chromosome coordinate indexing system.
psi$chr <-  sub("^M$","MT",sub("chr","",psi$chr))
psi$start <- psi$end <- psi$coor

psiz <- cbind(psi[,c("chr","start","end")],
              strand=".", info="",
              gene=psi$Annotation.first,
              psi=psi[,"mm.DirectMINUSmm.IVT"]/100, # % -> fraction
              source=paste0("fanari24pre","_",psi$cell))
PSI <- rbind(PSI, psiz)


### PROCESS PSI DATA - REMAP TO ENSEMBL MANE

## get ensembl protein ID
PSI$ensembl <- names(ens2nam[match(PSI$gene,  ens2nam)]) 

## take mean position - only affects @Zhang which provides ranges,
## and we found only one overlap (2nd and 3rd position of one AAS)
PSI$coor <- round(apply(PSI[,c("start","end")],1,mean))

## chromosome index
PSI$chr <-  chrIdx[PSI$chr]


### MAP PSI and AAS SITES
psidx <- coor2index(PSI[,c("chr","coor")], chrS=chrS)[,2]

## REANNOTATE PSI, EACH TO ENSEMBL MANE TRANSCRIPTS
cds.file <- "/home/raim/data/mammary/originalData/transcript_coordinates.tsv"
cds <- read.delim(cds.file)
cds <- cds[cds$transcript%in%genes$MANE,]
cdsi <- cds[,c("transcript","chr","start","end")] # SKIP STRAND!
cdsi$chr <- chrIdx[cdsi$chr]
cdsi <- coor2index(cdsi, chrS=chrS)

## find PSI sites in ensembl MANE transcript coordinates
## (ignoring strand)
## TODO: more efficient via cut or findInterval,see
## cut(x, breaks = interval.vector, include.lowest = TRUE)
psiu <- unique(psidx) # TODO: search only unique sites and map back
psi2ens <- sapply(psiu, function(x) which(cdsi$start <= x &  cdsi$end >= x))

psi2trans <- lapply(psi2ens,
                    function(x) {
                        if ( length(x)>0 ) {
                            x <- cdsi[x,"transcript"]
                        } else { x <- NA }
                        x
                    })


## expand back to site
PSI$transcripts <- unlist(lapply(psi2trans[match(psidx,psiu)],
                                function(x) paste(x, collapse=";")))

## take only first - NOTE: very few cases
PSI$transcript <- unlist(lapply(psi2trans[match(psidx,psiu)],
                                function(x) x[1]))
PSI$transcript[PSI$transcript=="NA"] <- NA


### ONLY CONSIDER PSI SITES THAT WE COULD MAP TO
### MANE ENSEMBL TRANSCRIPTS
cat(paste("NOTE: REMOVING", sum(is.na(PSI$transcript)), "of", nrow(PSI),
          "PSI SITES w/o MANE TRANSCRIPT MATCH\n"))

PSI <- PSI[!is.na(PSI$transcript),]
### RE-MAP PSI SITES COORDINATES
psidx <- coor2index(PSI[,c("chr","coor")], chrS=chrS)[,2]

## collect overlapping sites
psite <- usite
psia <- PSI
for ( i in 1:3 ) {

    bdidx <- get(paste0("bd",i,"dx"))
                 
    ## find same positions
    bid <- which(bdidx%in%intersect(psidx, bdidx))
    ##pid <- which(psidx%in%intersect(psidx, bdidx))
    cat(paste("found", length(bid), "overlapping sites", i, "\n"))
    
    ## annotate AAS
    psi2site <- match(bdidx, psidx)
    sitepsi <- PSI[psi2site, ]
    sitepsi$chr <-
        gsub("chrNA","",paste0("chr", chrMap[sitepsi$chr,2]))
    sitepsi$chr <- gsub(";","",sitepsi$chr)
    sitepsi$chr[sitepsi$chr==""] <- NA
    
    colnames(sitepsi) <- paste0("codon",i, "_", colnames(sitepsi))
    psite <- cbind(psite, sitepsi)

    ## annotate psi sites
    aas2psi <- match(psidx, bdidx)
    aasite <- psite[aas2psi,c("transcript","fromto","codon","RAAS.median")]
    colnames(aasite) <- paste0("codon",i, "_", colnames(aasite))
    psia <- cbind(psia, aasite)
}

## collapse all psi
## NOTE: using median as a helper

pcodons <- function(x) {
    y <- which(!is.na(x))
    if (length(y)==0) y<-NA;
    paste(y[!is.na(y)],collapse=";")
}
pvals <- function(x) paste0(x[!is.na(x)],collapse=";")
pmax <- function(x) {
    y <- x[!is.na(x)]
    ifelse(length(y)==0, NA, max(y))
}

## add summary to AAS table
psite$psi.n <- apply(psite[,paste0("codon",1:3,"_psi")], 1,
                     function(x) sum(!is.na(x)))
psite$psi.max <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pmax)
psite$psi <- apply(psite[,paste0("codon",1:3,"_psi")], 1, pvals)
psite$psi.codonpos <- apply(psite[,paste0("codon",1:3,"_psi")], 1,  pcodons)
psite$psi.source <- apply(psite[,paste0("codon",1:3,"_source")], 1, pvals)
psite$psi.chr <- apply(psite[,paste0("codon",1:3,"_chr")], 1, pvals)
psite$psi.coor <- apply(psite[,paste0("codon",1:3,"_coor")], 1, pvals)

## add summary to PSI table
psia$aas.n <- apply(psia[,paste0("codon",1:3,"_RAAS.median")], 1,
                    function(x) sum(!is.na(x)))
psia$aas.transcript <- apply(psia[,paste0("codon",1:3,"_transcript")], 1, pvals)



## CORRELATION OF PSI % TO RAAS

sources <- unique(na.omit(unlist(psite[,paste0("codon",1:3,"_source")])))
pch.source <- 1:length(sources) 
names(pch.source) <- sources

## for each source/cell line separately
for ( src in sources ) {
    ##ssite <- psite[grep(src,psite$psi.source),]
    src <- sub(".*_","",src)
    for ( i in 1:3 ) {

        
        ccol <- paste0("codon",i,"_psi")
        scol <- paste0("codon",i,"_source")

        ssite <- psite[grep(src,psite[,scol]),]
        nvals <- sum(!is.na(ssite[,ccol]))

        pc <- NA

        ## PLOT ONLY SIGNIFICANTLY CORRELATING
        if (nvals<3) next
        if ( cor.test(ssite$RAAS.median, ssite[,ccol])$p.value > .05 )
            next
        
        plotdev(file.path(rfig.path,paste0("psi_raas_codon",i,"_",src)),
                type=ftyp, height=2.5, width=2.5, res=200)
        par(mai=c(.5,.5,.25,.25), mgp=c(1.3,.3,0), tcl=-.25)
        if (nvals==0) {
            plot(1, col=NA, xlab=NA, ylab=NA, axes=FALSE)
            box()
        } else if ( nvals>2) {
            pc <- plotCor(ssite$RAAS.median, ssite[,ccol], density=FALSE,
                          xlab=NA, ylab=NA, signif=2, col=NA,
                          legcex=.8)$p
            
        } else {
            plot(ssite$RAAS.median, ssite[,ccol], col=NA,
                 xlab=NA, ylab=NA)
        }
        if ( nvals>0 ) {
            points(ssite$RAAS.median, ssite[,ccol],
                   pch=1,#pch.source[ssite[,scol]],
                   cex=.75, col=i)
            basicPlotteR::addTextLabels(ssite$RAAS.median, ssite[,ccol],
                                        labels=ssite$fromto, #xpd=TRUE,
                                        col.label=1, col.line=i,
                                        cex.label=.5)#, pos=4)
        }
        mtext(paste0(src, ", codon pos. ", i) , 3, 0)
        mtext(xl.raas, 1, 1.3)
        mtext(bquote(psi~fraction), 2, 1.3)
        dev.off()
    }
}

## split multiple psi per AAS
psil <- lapply(strsplit(psite$psi, ";"), as.numeric)
xy <- cbind(rep(psite$RAAS.median, lengths(psil)), unlist(psil))

plotdev(file.path(rfig.path,"psi_raas_codons"), type=ftyp,
        height=2.5, width=2.5, res=200)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(xy[,1], xy[,2], density=FALSE,
        xlab=xl.raas, ylab=expression(psi~fraction), pch=NA)
points(psite$RAAS.median, psite$codon1_psi, col=1, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon2_psi, col=2, pch=19, cex=.5)
points(psite$RAAS.median, psite$codon3_psi, col=3, pch=19, cex=.5)
legend("right",paste(1:3), title="codon position", col=1:3, cex=.5,
       pch=19, pt.cex=.75, bty="n")
dev.off()


### HYPERGEO, USING NUMBER OF Us AS BACKGROUND

## COUNT BACKGROUND: Us in CODING SEQUENCE
tfas <- readFASTA(tfas.file, grepID=TRUE)

pids <- unique(c(psite$transcript, PSI$transcript))
sfas <- tfas[pids[pids%in%names(tfas)]]

## total balls: all Us in coding sequences of the total gene set
totaa <- sum(unlist(lapply(sfas,
                           function(x) sum(unlist(strsplit(x$seq,""))=="T"))))

## white balls: all UNIQUE detected psi
m <- length(unique(paste(PSI$chr, PSI$coor)))
n <- totaa-m # black balls

## balls drawn: all U in AAS codons
## take only unique site here!
FILT <- !duplicated(paste(psite$chr,
                          psite$coor, psite$strand))
k <- sum(unlist(strsplit(psite$codon[FILT],""))=="T")

## q: white balls drawn number of AAS with psi
q <- sum(psite$psi.n[FILT])

##m: the number of white balls in the urn.
##n: the number of black balls in the urn.
##k: the number of balls drawn from the urn

p <- phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE)

## pdist

plotdev(file.path(rfig.path,"psi_hypergeotest"), res=200,
        width=2.5, height=2.5, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
qs <- 1:(2*q)
plot(qs, phyper(q=qs-1, m=m, n=n, k=k, lower.tail=FALSE),
     xlab="overlap count", ylab="probability", type="l", lwd=2)
legend("topright", c(paste("transcripts=", round(length(pids)/1000,1),"k"),
                     paste("total U=", round(totaa/1000000,1),"M"),
                     paste("psi sites=", round(m/1000,1),"k"),
                     paste("U in AAS=", round(k/1000,1),"k"),
                     paste("overlap=", q)), cex=.8,
       box.col=NA, bg=NA)
points(q, p, pch=4, cex=1.2, col=2)
shadowtext(q, p,label=paste0("p=",signif(p,1)), pos=3, col=2, font=2)
##figlabel("union",pos="bottomleft", font=2)
dev.off()

