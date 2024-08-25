
## GO ANALYSIS of PROTEINS WITH AMINO ACID SUBSTITUTIONS

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")


## only use unique BP/SAAP per Dataset
dfig.path <- file.path(fig.path,"function")
dir.create(dfig.path, showWarnings=FALSE)

tmtu <- tmtf
value <- "RAAS"


## TIGHT PLOT FOR MAIN
p.tgt <- 1e-15 # tighter cutoff for GO analysis

## BP/SAAP 
tmtu$BP.SAAP <- paste0(tmtu$BP,".",tmtu$SAAP)


### RAAS dotplot for GOslim CATEGORIES

## get GOslim table for use with clusterAnnotation
got <- parseAnnotationList(genes[,c("ID","GOslim")]) 
## replace GO IDs by terms
terms <- read.delim(goslim.file)
trms <- terms[,2]
names(trms) <- terms[,1]
colnames(got) <- trms[colnames(got)]

mid <- "ENSG00000277856"

## shorten names
colnames(got) <- sub("plasma membrane", "PM",
                     sub("binding", "bnd.",
                         sub("templated", "templ.",
                             sub("regulation", "reg.",
                                 sub("transcription", "transcr.",
                                     sub("localization","local.",
                                         colnames(got)))))))

## calculate RAAS profile
go.ovl <- raasProfile(x=tmtu, id="BP.SAAP", 
                    rows=got[tmtu$gene,], cols="Dataset",
                    bg=TRUE, value=value, 
                    col.srt=uds,
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=0)
## sort and filter significant
##ovc <- sortOverlaps(go.ovl, axis=2, p.min=p.min, sign=1, cut=TRUE)

## NOTE: remnant from former loop over several annotations
cids <- c("go")
names(cids) <- cids
cnms <- cids
cnms["go"] <- "GOslim"

omai <- c(.8,CMAIL,.6,.6)
for ( cid in cids ) {

    ovf <- get(paste0(cid, ".ovl"))
    did <- paste0("type_", cid, "_")

    ## p.tgt, high RAAS
    ovc <- sortOverlaps(ovf, p.min=p.tgt, cut=TRUE, sign=1)
    if ( nrow(ovc$p.value)>0 ) {        
        fname <- file.path(dfig.path,paste0(did,SETID,"_ptgt_high"))
        
        lmai <- omai
        lmai[2] <- .1 + .1*max(nchar(rownames(ovc$p.value)))
        plotProfiles(ovc, fname=fname,
                     mai=lmai, ttcols=ttcols, value="median",
                     dot.sze=dot.sze, p.dot=p.dot,
                     p.min=p.min, p.txt=p.txt,
                     rlab=ovc$p.min, llab=cnms[cid], ftyp=ftyp,
                     mtxt="", mtxt.line=2.3,
                     vcols=acols, vbrks=abrks,
                     gcols=gcols, ffam=FONT, bg=NA, plot.all=FALSE)
    }  
}

### EXTRACT HIGH RAAS PROTEINS FOR EACH GO

ovc <- sortOverlaps(go.ovl, p.min=p.tgt, cut=TRUE, sign=1)
## .. and MEASURED in all Datasets!
all <- rownames(ovc$count)[apply(ovc$count, 1, function(x) all(x>0))]
ovc <- sortOverlaps(ovc, srt=all, cut=TRUE)
ovc <- sortOverlaps(ovc, p.min=p.txt, cut=FALSE)#, sign=1)

## count unique proteins for each enriched GO
gol <- list()
for ( i in 1:nrow(ovc$p.value) ) 
    gol[[i]] <- rownames(got)[which(got[,rownames(ovc$p.value)[i]])]
## filter those in TMT set
gol <- lapply(gol, function(x) unique(x[x%in%tmtu$gene]))
names(gol) <- rownames(ovc$p.value)
## GET PROTEIN NAMES
gon <- lapply(gol, function(x) tmtu$name[tmtu$gene%in%x])


## PLOT GENES FROM GO ENRICHMENT CATEGORIES 

## calculate raas profiles for each gene (via name column)
pr.ovl <- raasProfile(x=tmtu, id="SAAP", 
                   rows="name", cols="Dataset",
                   col.srt=uds,
                   bg=TRUE, value=value, 
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)
for ( i in 1:length(gon) ) {

    goid <- names(gon)[i]

    ## get list of proteins that contribute to GO class
    gopr <- sortOverlaps(pr.ovl, axis=2, srt=unique(unlist(gon[[i]])))
    ## filter with a mild p-value cutoff, at high RAAS
    gopr <- sortOverlaps(gopr, p.min=p.txt, sign=1, cut=TRUE)
    ## sort and cut by number of RAAS values,
    nsrt <- names(sort(gopr$num.query[gopr$num.query[,1]>1,1]))
    gopr <- sortOverlaps(pr.ovl, axis=2, srt=nsrt)
    
    fname <- file.path(dfig.path,paste0("goslim_",SETID,"_",
                                        gsub(" ","_",goid)))
    omai <- c(.2,1.1,.05,.5)
    nw <- ncol(gopr$p.value)*.2 + omai[2] + omai[4]
    nh <- nrow(gopr$p.value)*.2 + omai[1] + omai[3]

    plotdev(fname,  height=nh, width=nw, res=300, type=ftyp)
    par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family="sans")
    dotprofile(gopr, value="median",
               vcols=acols, vbrks=abrks, p.dot=p.dot, dot.sze=dot.sze,
               axis=2,  show.total=TRUE, xlab=NA, ylab=NA)
    figlabel(goid, pos="bottomleft", font=2, cex=1)
    dev.off()
}
