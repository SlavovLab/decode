
## CODON FREQUENCY ANALYSIS of AMINO ACID SUBSTITUTION SITES

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")


## local output path
cfig.path <- file.path(fig.path,"codons")
dir.create(cfig.path, showWarnings=FALSE)

### START ANALYSIS


### CODON FREQUENCIES vs. RAAS

## ONLY USE RAAS WITH ASSIGNED CODONS
ctmt <- tmtf[tmtf$codon!="",]

## add columns for codon positions
cpos <- strsplit(ctmt$codon,"")
ctmt$pos1 <- unlist(lapply(cpos, function(x) x[1]))
ctmt$pos2 <- unlist(lapply(cpos, function(x) x[2]))
ctmt$pos3 <- unlist(lapply(cpos, function(x) x[3]))

## LOAD GLOBAL GENE-WISE CODON COUNTS
codons <- read.delim(codon.file, row.names=1)

## Wu et al. 2019: Codon Stability Coefficient
csc <- read.csv(wu19.file, row.names=1)


## CODON COUNT and FREQUENCY IN ALL UNIQUE TRANSCRIPTS 

## codon counts in all mapped transcripts
## NOTE: this is used as our background frequency
cod  <- codons[unique(ctmt$transcript),]
codt <- apply(cod,2,sum) # total count

if ( any(!ctmt$transcript%in%rownames(codons)) )
    stop("some codon frequencies not found, rerun calculation")

### CODON SORTING

## PER AA
## transcript codon frequencies per AA and SORTING by frequency and AA prop
codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
codl <- lapply(codl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
codl <- codl[aap.srt[aap.srt%in%names(codl)]] ## SORT BY AA PROP

## LOCAL CODON SORTING by background frequencies
## (AA property ->codon frequency)
codon.srt <- sub("\\.","-",names(unlist(codl)))
## remove K|R codons
codon.srt <- codon.srt[grep("[KR]", codon.srt, invert=TRUE)]

## codon frequencies per AA
Fbg <- lapply(codl, function(x) x/sum(x)) # codon frequency
Fbg <- unlist(Fbg)
names(Fbg) <- sub("\\.","-",names(Fbg))



## median codon RAAS:
## NOTE: median over all measurements
Craas <- sapply(COD.SRT, function(cl)
    log10(median(10^ctmt$RAAS[ctmt$aacodon==cl])))

## calculate unique SAAP/BP medians before
## taking codon median RAAS
## NOTE: THIS OVERRIDES the site table from init

csite <- split(ctmt$RAAS, ctmt$unique.site)
csite <- listProfile(csite, y=ctmt$RAAS, use.test=use.test, min=3)
site.codons <- lapply(split(ctmt$aacodon,  ctmt$unique.site), unique)
if ( length(table(lengths(site.codons)))!=1 )
    stop("multiple codons per protein site")
csite <- cbind(csite, codon=unlist(site.codons))

## median of site medians
## NOTE: correlation to codon frequency is lost
## when using the median of medians.
Craas.site <- sapply(COD.SRT, function(cl)
    log10(median(10^csite$median[csite$codon==cl])))

## add columns for codon positions
cpos <- strsplit(sub("[A-Z]-","",csite$codon),"")
csite$pos1 <- unlist(lapply(cpos, function(x) x[1]))
csite$pos2 <- unlist(lapply(cpos, function(x) x[2]))
csite$pos3 <- unlist(lapply(cpos, function(x) x[3]))


## CODON FREQUENCIES at AAS
lcodt <- table(csite$codon)
names(lcodt) <- sub(".*-","",names(lcodt))
## AAS codon frequencies per AA
lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])

lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
lcodl <- lcodl[names(aaprop)] ## ##SORT BY AA PROP

## store for sanity check:
## why more codons via raasProfile unique counter than here?
lcodt.global <- lcodt
Fbg.lst <- codl
Faas.lst <- lcodl

## AA-specific codon frequency
Faas <- unlist(lapply(lcodl, function(x) x/sum(x)))
names(Faas) <- sub("\\.","-", names(Faas))

## inspect
if ( interactive() ) { # properly plotted below loop
    plotCor(Craas.site[names(Craas)], Craas)
    plotCor(Craas.site[names(Faas)], Faas)
    plotCor(Craas[names(Faas)], Faas)
    plotCor(Craas[names(Fbg)], Fbg)
    plotCor(Fbg[names(Faas)], Faas, xlim=c(0,1), ylim=c(0,1))
}

## codon rank - TODO: remove or improve this
codr <- lapply(codl, function(x) round(6*rank(x)/length(x))) # codon rank
codr$W[1] <- 0
codr$M[1] <- 0
## codon rank
codr <- unlist(codr)
names(codr) <- sub("\\.", "-", names(codr))
ctmt$codon.rank <- as.character(codr[ctmt$aacodon])
    
## codon frequency bins
ctmt$codon.bins <- cut(Fbg[ctmt$aacodon], seq(0,1,.1), include.lowest = TRUE)

## TODO: do by 2|3, 4 or 6 codons
cod.bins <- cut(Fbg[ctmt$aacodon], seq(0,1,.1), include.lowest = TRUE)
ncod <- unlist(lapply(CODL, length))


## PER CODON RAAS PROFILES per Dataset

## matrix to store codon frequencies in bg and at AAS
allcodons <- sub("\\.","-",names(unlist(codl)))
codon.results <- matrix(NA, nrow=length(COD.SRT), ncol=length(auds))
colnames(codon.results) <- auds
rownames(codon.results) <- COD.SRT
codon.Fbg <- codon.Faas <- codon.raas <- codon.results


## CODON PLOT SUPPLEMENT - dotplots and frequencies
ds <- "all"

dtmt <- ctmt
dsLAB <- ""

## RAAS PROFILE ALL DATA
## NOTE: bg=FALSE, no column-wise background !
dtmt$all <- "all"
ova <- raasProfile(x=dtmt, id="unique.site", 
                   rows="all", cols="aacodon",
                   col.srt=codon.srt, filter=FALSE,
                   bg=FALSE, use.test=use.test, do.plots=FALSE, 
                       verb=0)

## RAAS PROFILE by Dataset
## NOTE: bg=TRUE :column-wise background
ovd <- raasProfile(x=dtmt, id="unique.site", 
                       rows="Dataset", cols="aacodon",
                       col.srt=codon.srt, filter=FALSE,
                       bg=TRUE, bg.dir="row",
                       use.test=use.test, do.plots=FALSE, 
                   verb=0)

## sort rows
nr <- uds[uds%in%rownames(ovd$p.value)]
if ( length(nr)>1 )
    ovd <- sortOverlaps(ovd, axis=2, srt=nr)



## RAAS PROFILE by Incorporated AA
## NOTE: bg=FALSE, no column-wise background !
ovw <- raasProfile(x=dtmt, id="unique.site", 
                   rows="to", cols="aacodon",
                   col.srt=codon.srt, filter=FALSE,
                   bg=FALSE, use.test=use.test, do.plots=FALSE, 
                   verb=0)

## sort by amino acid property via shapely colors!
aasrt <- aap.srt[aap.srt%in%rownames(ovw$p.value)]
ovw <- sortOverlaps(ovw, axis=2, srt=aasrt, cut=TRUE)


## CODON FREQUENCIES

dsite <- split(dtmt$RAAS, dtmt$unique.site)
dsite <- listProfile(dsite, y=dtmt$RAAS, use.test=use.test, min=3)
dsite.codons <- lapply(split(dtmt$aacodon,  dtmt$unique.site), unique)
if ( length(table(lengths(dsite.codons)))!=1 )
    stop("multiple codons per protein site")
dsite <- cbind(dsite, codon=unlist(dsite.codons))

## codon frequencies at AAS
lcodt <- table(dsite$codon)
names(lcodt) <- sub(".*-","",names(lcodt))
  
## AAS codon frequencies per AA
lcodl <- split(lcodt, GENETIC_CODE[names(lcodt)])

lcodl <- lapply(lcodl, sort, decreasing=TRUE) ## SORT BY MOST FREQUENT
lcodl <- lcodl[aap.srt] ## ##SORT BY AA PROP

## AA-specific codon frequency
Faas.ds <- unlist(lapply(lcodl, function(x) x/sum(x)))
names(Faas.ds) <- sub("\\.","-", names(Faas.ds))

## CODON COUNTS IN ALL UNIQUE TRANSCRIPTS
## DATASET-SPECIFIC
cod  <- codons[unique(dtmt$transcript),]
codt <- apply(cod, 2, sum)

## transcript codon frequencies per AA
codl <- split(codt, GENETIC_CODE[sub(".*\\.","",names(codt))])
codl <- lapply(codl, sort, decreasing=TRUE) # SORT BY MOST FREQUENT
codl <- codl[aap.srt[aap.srt%in%names(codl)]] # SORT BY AA PROP
Fbg.ds <- unlist(lapply(codl, function(x) x/sum(x)))
names(Fbg.ds) <- sub("\\.","-", names(Fbg.ds))

## CODON FREQUENCIES
Fbg.ds.lst <- codl
Faas.ds.lst <- lcodl

## MEDIAN CODON RAAS
## get median RAAS for each codon and plot against codon frequency
Craas.ds <- sapply(codon.srt, function(cl)
    log10(median(10^dtmt$RAAS[dtmt$aacodon==cl])))

codon.Fbg[names(Fbg.ds),ds] <- Fbg.ds
codon.Faas[names(Faas.ds),ds] <- Faas.ds
codon.raas[names(Craas.ds), ds] <- Craas.ds

    
## CODON SORTING and PLOT SETTINGS

## sort by background (AA property ->codon frequency)
## i.e. as codl
cdsrt <- sub("\\.","-",names(unlist(codl)))

## FILTER AVAILABLE
cdsrt <- cdsrt[cdsrt%in%colnames(ovw$p.value)]

## codon frequency bins
codon.bins <- cut(Faas.ds[dtmt$aacodon], seq(0,1,.2), include.lowest = TRUE)

## codon class lines in plots
cdn.cnt <- table(sub("-.*","",colnames(ovw$p.value)))[aap.srt]
cdn.cnt <- cumsum(cdn.cnt[!is.na(cdn.cnt)])
    
## SORT CODON FREQUENCIES as in main plots
## TODO: constant sorting of  ALL codons 
Faas.ds <- Faas.ds[cdsrt]
Fbg.ds <- Fbg.ds[cdsrt]
Craas.ds <- Craas.ds[cdsrt]


plotProfiles(ova, fname=file.path(cfig.path,paste0("codon_",
                                                   SETID,"_",ds,
                                                   "_all")),
             fw=.2, mai=c(.05,.75,.05,.6), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, llab=dsLAB, ftyp=ftyp, mtxt="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, col.lines=cdn.cnt, ffam="monospace")


plotProfiles(ovd, fname=file.path(cfig.path,paste0("codon_",
                                                   SETID,"_",ds,"_Dataset")),
             fw=.2, mai=c(.05,.75,.05,.6), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, llab=dsLAB, ftyp=ftyp, mtxt="",
             vcols=vcols, vbrks=vbrks,
             gcols=gcols, col.lines=cdn.cnt, ffam="monospace")
maimain <- c(.6,.75,.6,.6) # re-used for codon axis
plotProfiles(ovw, fname=file.path(cfig.path,paste0("codon_",SETID,"_",ds)),
             fw=.2, mai=maimain, ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB, llab=dsLAB, ftyp=ftyp,
             mtxt="Incorporated AA", mtxt.cex=1.5,
             vcols=vcols, vbrks=vbrks,
             axis1.col=aap.cols[sub("-.*","",colnames(ovw$p.value))],
             axis2.col=aap.cols[rownames(ovw$p.value)],
             gcols=gcols, col.lines=cdn.cnt, ffam="monospace", verb=1)

### CODON FREQUENCY ANALYSIS


## full codon dotplots!
mai <- c(.05,.75,.1,.6)
fw <- .2
nw <- ncol(ovw$p.value) *fw + mai[2] + mai[4]


## codon plot, bottom row
plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",ds,
                                   "_codons_frequency_raas")),
        type=ftyp, res=300, width=nw,height=2.9)
par(mai=c(1,.75,.05,.6), mgp=c(1.3,.3,0), tcl=-.25, xaxs="i",
    family="monospace")
tmp <- boxplot(dtmt$RAAS ~ factor(dtmt$aacodon, levels=cdsrt),
               las=2, xlab="", ylab=xl.raas, axes=FALSE,
               cex=.5, pch=19, pars=list(outcol="#00000055"))
##axis(3, at=1:length(cdsrt), labels=round(Faas.ds[cdsrt],2), las=2)
axis(4); axis(2)
Map(axis, 1, at=1:length(tmp$names), labels=tmp$names, las=2,
    col.axis=aap.cols[sub("-.*","",tmp$names)], col=NA, font=2)
abline(v=.5+cdn.cnt, lwd=1, xpd=FALSE)
axis(1, at=c(0,.5+cdn.cnt), labels=FALSE, tcl=-.5)
figlabel(dsLAB, pos="bottomleft", font=2, cex=1.2)
figlabel(LAB, pos="bottomright", cex=.7)
text(c(3.25,13,25,42),rep(-9.5, 5), cex=1.2, font=2,
     labels=aaprop.srt, col=aaprop.cols[aaprop.srt], xpd=TRUE)
mtext("mRNA codons at amino acid substitution sites", 1, 4, cex=1.5)
dev.off()



## codon plot, row 1
mai.bar <- mai
mai.bar[c(2,4)] <- mai.bar[c(2,4)] +.05
plotdev(file.path(cfig.path,paste0("codon_",SETID,"_",
                                   ds,"_codons_barplot_onlyfbg")),
        type=ftyp, res=300, width=nw, height=1)
par(mai=mai.bar, mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp  <- barplot(rbind(Fbg.ds[names(Faas.ds)]),
               axes=FALSE, xlab=NA, beside=TRUE,
               ylab="", las=2, xaxt='n')
cdn.cn <- head(cdn.cnt, length(cdn.cnt)-1)
abline(v=bp[,cdn.cn]+unique(diff(bp[1,]))/2)

for ( ax in c(2,4) ) {
    axis(ax, at=c(0,.5,1), las=2, line=.2)
    axis(ax, at=seq(0,1,.1), labels=FALSE, line=.2, tcl=par("tcl")/2)
}
mtext("codon\nfrequency", 2, 2)
dev.off()

## CODON FREQUENCY ANALYSIS

## RESULT TABLE
codon.freq <- data.frame(codon=sub(".*-","",names(Fbg)),
                         AA=sub("-.*","",names(Fbg)),
                         frequency_all=Fbg,
                         frequency_AAS=Faas[names(Fbg)],
                         CSC=csc[sub(".*-","",names(Fbg)),"X293T_endo"],
                         RAAS=Craas[names(Fbg)])

write.table(codon.freq, file=file.path(out.path, "codon_frequencies.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, na="")

fftyp <- ftyp # "pdf" # 

## CODON RAAS vs. FREQUENCY
plotdev(file.path(cfig.path,paste0("codons_raas_fbg")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1), ##outliers=which(Fbg==1),
        line.methods=c("ols"), density=FALSE, pch=19,
        col="#00000099",
        xlab=expression("Relative codon frequency"~f[bg]), ylab=NA)
mtext(xl.raas, 2, 1.2)
dev.off()

## CODON RAAS vs. CSC
## Wu et al. 2019: Codon Stability Coefficient
plotdev(file.path(cfig.path,paste0("codons_raas_wu19")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.4,.3,0), tcl=-.25)
plotCor(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
        line.methods=c("ols"), density=FALSE, ylab=NA,
        xlab=NA, col=NA)
mtext("Codon Stability Coefficient", 1, 1.2)
points(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
       lwd=NA, cex=1, pch=19, col="#00000099")
mtext(xl.raas, 2, 1.2)
dev.off()


## encoded amino acid for coloring
aatmp <- sub("-.*", "", names(Fbg))
aasort <- names(sort(aa.cols))
aasort <- aasort[aasort%in%aatmp]

plotdev(file.path(cfig.path,paste0("codons_raas_fbg_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(Fbg, Craas[names(Fbg)], xlim=c(0,1),
        line.methods=c("tls","ols"), density=FALSE, pch=19,
        xlab=expression("Relative codon frequency"~f[bg]),
        ylab=xl.raas, col=NA)
points(Fbg, Craas[names(Fbg)], col=aa.cols[aatmp], pch=aa.pchs[aatmp], lwd=2)
dev.off()

plotdev(file.path(cfig.path,paste0("codons_raas_wu19_byAA")),
        type=fftyp, res=300, width=3,height=3)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
        line.methods=c("ols","tls"), density=FALSE, ylab=xl.raas,
        xlab="Codon Stability Coefficient", col=NA)
points(csc[sub(".*-","",names(Craas)),"X293T_endo"], Craas, 
       lwd=2, cex=1, col=aa.cols[sub("-.*","",names(Craas))],
       pch=aa.pchs[sub("-.*","",names(Craas))])
dev.off()


