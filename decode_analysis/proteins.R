
## PROTEIN-LEVEL ANALYSIS OF AMINO ACID SUBSTITUTIONS

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")

## local output path
pfig.path <- file.path(fig.path,"proteins")
dir.create(pfig.path, showWarnings=FALSE)


## figure settings
corW <- corH <- 2.5
pmai <- c(.5,.5,.25,.25)
pmpg <- c(1.3,.3,0)

## maximal SAAP/peptide, etc.
max.saap <- 12

## axis labels
xl.hlfm <- expression(protein~"half-life"/h)
xl.hlf <- expression(protein~"half-life"/h)

## overrule specific y-axis label
xl.prota <- xl.raas



### START ANALYSIS

### BP/SAAP RAAS and STATS

## bp median raas
bpl <- split(tmtf$RAAS, tmtf$BP)
bpstat <- listProfile(bpl, y=tmtf$RAAS, use.test=use.test, min=3)
## add data
bpstat$bplen <- nchar(rownames(bpstat))
## count distinct SAAP per BP
bps <- unlist(lengths(lapply(split(bdat$SAAP, bdat$BP), unique)))
bpstat$nsaap <- bps[rownames(bpstat)]

## add median intensities per BP+SAAP
bsints <- tmtf$BP.abundance+tmtf$SAAP.abundance
bpint <- unlist(lapply(split(bsints, tmtf$BP), median, na.rm=TRUE))

bpstat$intensity <- bpint[rownames(bpstat)]


### PROTEIN RAAS AND STATS

## protein median raas per protein w/o site-specific median first
ptl <- split(tmtf$RAAS, tmtf$ensembl)
ptstat <- listProfile(ptl, y=tmtf$RAAS, use.test=use.test, min=3)

## count distinct BP per  protein
pts <- unlist(lengths(lapply(split(bdat$BP, bdat$ensembl), unique)))
ptstat$nbp <- pts[rownames(ptstat)]
## count distinct SAAP per protein
pts <- unlist(lengths(lapply(split(bdat$SAAP, bdat$ensembl), unique)))
ptstat$nsaap <- pts[rownames(ptstat)]

## rank proteins by RAAS
ptstat$rank <- rank(ptstat$median)

### COLLECT PROTEIN DATA

#### TODO: why do we have n=747 at exon but n=713 for
#### the log intensity plot on exon/intron?

## PROTEIN ABUNDANCE
ptstat$intensity <- pint[rownames(ptstat)]

## PROTEIN LENGTH
ptstat$length <- plen[rownames(ptstat)]


## PROTEIN HALF-LIFE @Mathieson2018
ptstat$halflife <- phlv[match(pnms[rownames(ptstat)], names(phlv))]


## GET WHOLE PROTEIN MEAN IUPRED3 SCORE (from saap_peptides4.tsv
iu3 <- split(hdat$iupred3.protein, hdat$ensembl)
## QC: all protein level
table(unlist(lapply(iu3, function(x) length(unique))))
iu3 <- unlist(lapply(iu3, unique))

ptstat$iupred3 <- iu3[rownames(ptstat)]



## write-out collected protein results
write.table(cbind(protein=rownames(ptstat), ptstat),
            file=file.path(out.path, "proteins_raas.tsv"), sep="\t",
            na="", row.names=FALSE, quote=FALSE)


### AAS CLUSTERS: RAAS DISTRIBUTIONS of SINGLE BP
pns <- bpstat$nsaap
pns[pns>max.saap] <- paste0(">",max.saap)
paas <- table(pns)[c(as.character(1:max.saap),paste0(">", max.saap))]

bid <- rownames(bpstat)[which(bpstat$nsaap>80)]

plotdev(file.path(pfig.path,paste0("hotspots_SAAP_per_peptide_log")),
            type=ftyp, width=3, height=3, res=200)
par(mai=c(0.5,.5,.15,.05),mgp=c(1.3,.3,0), tcl=-.25, xaxs="i")
bp <- barplot(paas, xlab="SAAPs / base peptide",ylab=NA,
              axes=FALSE, axisnames = FALSE, log="y", ylim=c(1,max(paas)))
axis(1, at=bp, labels=names(paas), las=2, cex=.8)
options(scipen=999)
axis(2, at=1)
axis(2, at=10^(1:10))#, las=1)
axis(2, at=c(1:10,(1:10)*10,(1:10)*100, (1:10)*1000), tcl=-.125,
     labels=FALSE)
options(scipen=0)
mtext("Number of base peptides", 2, 1.3)
legend("topright",
       paste(sum(bpstat$nsaap>1,nr.rm=TRUE), "base peptides\nwith >1 SAAP"),
       bty="n")
dev.off()




### COMPARE PROTEIN RAAS TO VARIOUS PROTEIN LEVEL MEASURES

### PROTEIN INTENSITIES

plotdev(file.path(pfig.path,paste0("protein_intensities_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$intensity), ptstat$median, ylab=xl.raas,
        xlab=expression(log[10](intensity)), title=TRUE, cor.legend=FALSE)
dev.off()


### PROTEIN LENGTH

plotdev(file.path(pfig.path,paste0("protein_lengths_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$length), ptstat$median,
        ylim=range(ptstat$median), axes=FALSE,
        ylab=xl.prota, xlab="protein length", legpos="bottomright",
        legbg="#ffffff77", title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()


## PROTEIN HALF-LIVES, @Mathieson2018

## halflives all
plotdev(file.path(pfig.path,paste0("protein_halflives_all")),
        type=ftyp, res=300, width=corW,height=corH)
par(mai=pmai, mgp=pmpg, tcl=-.25)
plotCor(log10(ptstat$halflife), ptstat$median, ylim=range(ptstat$median),
        ylab=xl.prota, xlab=xl.hlfm, axes=FALSE, title=TRUE, cor.legend=FALSE)
axis(2)
axis(1, at=1:10, labels=10^(1:10))
axis(1, at=log10(rep(1:10, 5) * 10^rep(0:4, each=10)), tcl=-.125, labels=FALSE)
##box()
dev.off()






