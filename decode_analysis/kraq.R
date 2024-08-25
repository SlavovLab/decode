
### ANALYZE SEQUENCE CONTEXT along BASE PEPTIDES

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")

## local output path
kfig.path <- file.path(fig.path,"kraq")
dir.create(kfig.path, showWarnings=FALSE)

### remove duplicate sites!
cdat <- bdat[!duplicated(paste(bdat$ensembl,bdat$pos, bdat$fromto)),]

### POSITION OF AAS IN PEPTIDES

## Categorize by N+i and C-i
mx <- 9

nt <- cdat$site
ct <- nchar(cdat$BP)-cdat$site +1

## fuse central AAS
nt[nt>mx] <- ct[ct>mx] <- mx+1

## take smaller
at <- paste0("N",nt)
at[ct<nt] <- paste0("C",ct[ct<nt])

## re-name central AAS
at[at==paste0("N",mx+1)] <- paste0(">",mx)
at.srt <- c(
    paste0("N",1:mx),
    paste0(">",mx),
    paste0("C",mx:1))

asite.bins <- at
asite.srt <- at.srt

## calculate overlap enrichment by cumulative hypergeometric
## distribution tests
cls <- clusterCluster(cdat$fromto, asite.bins, cl2.srt=asite.srt)

## sort by only enriched
cls <- sortOverlaps(cls, p.min=1e-3)
clc <- sortOverlaps(cls, p.min=p.txt, cut=TRUE, sign=1)

plotdev(file.path(kfig.path,paste0("peptides_AAS_AAStype_tight")),
        height=.2*nrow(clc$p.value)+1.1, width=.25*nrow(clc$p.value)+1.25,
        res=300, type=ftyp)
par(mai=c(.5,.75,.5,.5), mgp=c(2,.3,0), tcl=-.05, family="monospace")
plotOverlaps(clc, p.min=p.min, p.txt=p.txt, ylab=NA, xlab=NA,
             show.total=TRUE, show.sig=FALSE, axis=NA, text.cex=.7)
par(mgp=c(1.3,.3,0), tcl=-.25)
axis(1, at=1:ncol(clc$p.value), labels=colnames(clc$p.value), las=2)
axex <- ftlabels(rownames(clc$p.value))
axis(2, length(axex):1, labels=axex, las=2, family="monospace")
mtext("AAS type", 2, 2.7)
mtext("position of AAS in peptide", 1, 1.5)
##figlabel("AAS", pos="bottomright",cex=1.2, font=2)
dev.off()
