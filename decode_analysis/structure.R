
## BP/SAAP-LEVEL CORRELATION OF CONSERVATION AND DISORDER

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("tmtf") )
    source("raas_init.R")

## local output path
sfig.path <- file.path(fig.path,"structure")
dir.create(sfig.path, showWarnings=FALSE)

### START ANALYSIS


fname <- file.path(sfig.path,paste0("structure_cor_MMSeq2_RAAS"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$MMSeq2, bdat$RAAS,
        xlab="conservation, MMSeq2",
        ylab=xl.raas, legpos="topright")
dev.off()


fname <- file.path(sfig.path,paste0("structure_cor_iupred3_RAAS"))
plotdev(fname, height=3.5, width=3.5, res=300, type=ftyp)
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.3,0), tcl=-.25)
plotCor(bdat$iupred3, bdat$RAAS,
        xlab="disordered score, iupred3",
        ylab=xl.raas, legpos="bottomleft")
dev.off()


