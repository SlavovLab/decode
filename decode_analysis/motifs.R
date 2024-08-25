
## SEQUENCE CONTEXT at AMINO ACID SUBSTITUTION SITES

## project-specific functions
source("raas_utils.R")

## common initialization of BP/SAAP mapping and TMT level RAAS data
## loading, mapping, filtering, data selection, output paths,
## ID mappings, etc.
if ( !exists("bdat") )
    source("raas_init.R")

## local output path
mfig.path <- file.path(fig.path,"motifs")
dir.create(mfig.path, showWarnings=FALSE)

## AA SORTING
AAS <- sort(unique(GENETIC_CODE)) 
AAT <- AAS[AAS!="*"]
ABC <- DiffLogo::ASN 


## motif p-value cutoffs *, **, ***
psig <- 10^-c(3,5,10)

## remove identical sequences before motif plots
remove.duplicated.sequences <- TRUE

#### GENERAL SEQUENCE LOGOS

### remove duplicate sites!
## TODO: is there a cleaner way to do this? here we randomly loose AAS types
## TODO: is this still required? in logo calculation we rm duplicated anyways,
## but in th enrichment plots it may have effects
cdat <- bdat[!duplicated(paste(bdat$ensembl,bdat$pos, bdat$fromto)),]


### AA MATRIX
aam <- do.call(rbind, strsplit(cdat[["AA"]],""))
nc <- (ncol(aam)-1)/2
colnames(aam) <- -nc:nc
rownames(aam) <- paste0(cdat$BP,"_", cdat$SAAP)

## AA MATRIX WITH AAS
aas <- aam
aas[,"0"] <- cdat$to



## TODO: AAS_overlap for RAAS BINS

### HYPERGEO TESTS - tigther context

omai <- c(.5,.5,.6,.6)

ovl <- aaProfile(aam[,as.character(-7:7)], abc=AAT, alternative="greater")
ovl <- sortOverlaps(ovl, p.min=p.txt, sign=1)

nw <- ncol(ovl$p.value)*.25 + omai[2] + omai[4]
nh <- nrow(ovl$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("AAS_overlap")),
        height=nh, width=nw, res=300, type=ftyp)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
plotOverlaps(ovl, p.min=p.min, p.txt=p.txt, show.total=TRUE,
             xlab="Distance from AAS", ylab="Encoded AA",
             text.cex=.7, show.sig=FALSE)
##figlabel("all", pos="bottomleft")
dev.off()


### DEFINE MOTIF CLASSES
do.all.motifs <-  TRUE # FALSE # 

CTXT <- as.character(c(-2,-1,1,2))
CCxCC <- apply(aam[,CTXT], 1, function(x) any(x%in%c("C")))
MMxMM <- apply(aam[,CTXT], 1, function(x) any(x%in%c("M")))
WWxWW <- apply(aam[,CTXT], 1, function(x) any(x%in%c("W")))
GGxGG <- apply(aam[,CTXT], 1, function(x) any(x%in%c("G")))

classes <- cbind(
    "KRAQ"=cdat$to%in%c("G","A") & cdat$site%in%1:3,
    CCxCC=CCxCC,
    MMxMM=MMxMM,
    WWxWW=WWxWW,
    GGxGG=GGxGG
)

## NA values in filters
classes[is.na(classes)] <- FALSE

## store manual selection before (otpionally)
## expanding to full scan of AAS types
selected <- colnames(classes)

## ADD Encoded->Incorporated as CLASSES
fts <- unique(cdat$fromto)
ftcls <- matrix(NA, ncol=length(fts), nrow=nrow(aam))
colnames(ftcls) <-  paste0("fromto_",fts)
for ( i in seq_along(fts) ) 
    ftcls[,i] <- cdat$fromto==fts[i]
if ( !interactive() | do.all.motifs ) classes <- cbind(classes, ftcls)

frm <- unique(cdat$from)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("from_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- cdat$from==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

frm <- unique(cdat$to)
ftcls <- matrix(NA, ncol=length(frm), nrow=nrow(aam))
colnames(ftcls) <- paste0("to_",frm)
for ( i in seq_along(frm) ) 
    ftcls[,i] <- cdat$to==frm[i]
if ( !interactive() | do.all.motifs  ) classes <- cbind(classes, ftcls)

cols <- as.character(c(-3:3))

## custom ranges
rngs <- rep(list(cols), ncol(classes))
names(rngs) <- colnames(classes)


rngs$MMxMM <- rngs$WWxWW <- rngs$CCxCC <- rngs$CCxPP <-as.character(-2:2)
rngs$"KRAQ" <- as.character(-3:0)
## suppress p-value indicators
opval <- rep(list(""), ncol(classes))
names(opval) <- colnames(classes)
opval$MMxMM <- opval$WWxWW <- opval$CCxCC <- as.character(c(-2,-1,1,2))

## suppress on incorporated side only
npval <- opval
npval["KRAQ"] <- "0"

## use our internal AA colors
ABC$cols <- aa.cols[ABC$chars]
ABC$cols["V"] <- "#2eb774"


log.path <- file.path(mfig.path, "logos")
dir.create(log.path, showWarnings=FALSE)
tmp.path <- file.path(mfig.path, "selected")
dir.create(tmp.path, showWarnings=FALSE)

## PLOT ALL DIFFERENCE LOGOS
for ( i in 1:ncol(classes) ) {

    id <- colnames(classes)[i]
    lb <- id # sub(".*_", "", id)
    filt <- classes[,id]
    cols <- rngs[[id]]
    axlab <- as.character(cols)
    axlab[axlab=="0"] <- ""
 
    if ( sum(filt)<2 ) next

    tmp.path <- log.path
    if ( id %in% selected )
        tmp.path <- file.path(mfig.path, "selected")
    

    ## PWM ENCODED
    ## position weight matrices and diffLogo

    aam1 <- aam[ filt,,drop=FALSE]
    aam2 <- aam[!filt,,drop=FALSE]

    if ( remove.duplicated.sequences ) {
        aam1 <- aam1[!duplicated(apply(aam1,1,paste,collapse="")),]
        aam2 <- aam2[!duplicated(apply(aam2,1,paste,collapse="")),]
    }

    ## filter tight context
    aam1 <- aam1[,cols]
    aam2 <- aam2[,cols]

    pfm1 <- getPFM(aam1, alphabet=ABC$chars)
    pfm2 <- getPFM(aam2, alphabet=ABC$chars)

    n1 <- nrow(aam1)
    n2 <- nrow(aam2)
    
    dfob <- createDiffLogoObject(pfm1, pfm2, alphabet=ABC)
    dfop <- enrichDiffLogoObjectWithPvalues(dfob, n1=n1, n2=n2)

    ## suppress p.vals where ALL AA are equal
    ## -> obvious selection bias
    setpo <- apply(pfm1, 2, function(x) all(x%in%c(0.0,1.0)))
    ## manual pval suppression
    setpo[cols%in%opval[[id]]] <- TRUE
    dfop$pvals[setpo] <- 1

    ## PWM INCORPORATED
    ## position weight matrices and diffLogo
    aas1 <- aas[ filt,,drop=FALSE]
    aas2 <- aas[!filt,,drop=FALSE]
    
    if ( remove.duplicated.sequences ) {
        aas1 <- aas1[!duplicated(apply(aas1,1,paste,collapse="")),]
        aas2 <- aas2[!duplicated(apply(aas2,1,paste,collapse="")),]
    }
    
    ## filter tight context
    aas1 <- aas1[,cols]
    aas2 <- aas2[,cols]
    
    pfn1 <- getPFM(aas1, alphabet=ABC$chars)
    pfn2 <- getPFM(aas2, alphabet=ABC$chars)
    
    n1 <- nrow(aas1)
    n2 <- nrow(aas2)
    
    dfnb <- createDiffLogoObject(pfn1, pfn2, alphabet=ABC)
    dfnp <- enrichDiffLogoObjectWithPvalues(dfnb, n1=n1, n2=n2)
        
    ## suppress p.vals where ALL AA are equal
    ## -> obvious selection bias
    setpn <- apply(pfn1, 2, function(x) all(x%in%c(0.0,1.0)))
    ## manual pval suppression
    setpn[cols%in%npval[[id]]] <- TRUE
    dfnp$pvals[setpn] <- 1
    
       
    ## PLOT LOGOS

    if ( length(grep(":",lb)) ) {
        ft <- unlist(strsplit(lb,":"))
        ft[ft=="x"] <- ""
           lb <- as.expression(bquote(.(ft[1]) %->% .(ft[2])))
    } else lb <- as.expression(bquote(.(lb)))

    ## set ylim to non-prespecified positions (via setp/setn)
    focus.ylim <- FALSE
    if ( length(grep("^from", colnames(classes)[i]))>0 )
        focus.ylim <- TRUE

    ## ENCODED
    if ( any(dfop$pvals<min(psig)) | id %in% selected ) {

        mmai <- c(.5,.5,.25,.1)
        wd <- .4*length(cols) + mmai[2] + mmai[4]
        ht <- 2.5
        mn <- dfob$ylim.negMax
        mx <- dfob$ylim.posMax
        if ( focus.ylim )
            mx <- max(dfob$ymaxs[!setpo])
        
        
        plotdev(file.path(tmp.path,paste0("AA","_logos_", id,"_encoded")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=0, ymax=mx)
        axis(2)
        axis(1, at=1:length(cols), labels=axlab)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=2)
        figlabel(paste0("n=",nrow(aam1)), pos="topright", font=2)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT)
        diffLogo_addPvals(dfop, ymin=mx, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        dev.off()
    }


    ## INCORPORATED
    if ( any(dfnp$pvals<min(psig))  | id %in% selected ) {

        mmai <- c(.5,.5,.25,.1)
        wd <- .4*length(cols) + mmai[2] + mmai[4]
        ht <- 2.5
        mn <- dfnb$ylim.negMax
        mx <- dfnb$ylim.posMax
        if ( focus.ylim )
            mx <- max(dfnb$ymaxs[!setpn])
        
        plotdev(file.path(tmp.path,paste0("AA","_logos_", id,"_incorporated")),
                height=ht, width=wd, res=300, type=ftyp)
        par(mai=mmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfnb, sparse=TRUE, ymin=0, ymax=mx)
        axis(1, at=1:length(cols), labels=axlab)
        axis(2)
        mtext("JS divergence", 2, 1.3)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels=expression(INC), las=2)
        figlabel(paste0("n=",nrow(aas1)), pos="topright", font=2)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT)
        diffLogo_addPvals(dfnp, ymin=mx, levels=psig)
        ##mtext("incorporated", 4,-.25,adj=.05)
        dev.off()
    }
    
    ## TIGHT LOGO PLOT FOR MAIN
    if ( (any(dfop$pvals<min(psig))  |
          any(dfnp$pvals<min(psig))) |
         id %in% selected ) {
        
        mno <- dfob$ylim.negMax
        mxo <- dfob$ylim.posMax
        mnn <- dfnb$ylim.negMax
        mxn <- dfnb$ylim.posMax

        if ( focus.ylim ) {
            mxo <- max(dfob$ymaxs[!setpo])
            mxn <- max(dfnb$ymaxs[!setpn])
        }
        
        ##mxo <- mxn <- max(mxo, mxn)
        
        ## figure heights
        ht <- 3
        omai <- c(.15,.3,.25,.05)
        nmai <- c(.25,.3,.15,.05)

        ## heights
        htn <- ht/2 + nmai[1] + nmai[3]
        hto <- ht/2 + omai[1] + omai[3]
        
        ## adaptive width
        wd <- .4*length(cols) + nmai[2] + nmai[4]
        
        dfop$ylab <- "JS divergence"

        plotdev(file.path(tmp.path,paste0("AA","_logos_", id, "")),
                height=hto+htn , width=wd, res=300, type=ftyp)
        layout(t(t(1:2)), heights=c(hto, htn))

        ## encoded
        par(mai=omai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfob, sparse=TRUE, ymin=0, ymax=mxo)
        ## TODO: custom axis
        xat <- pretty(par("usr")[3:4])
        xlb <- c("0",xat[length(xat)-2])
        xat <- c(0,xat[length(xat)-2])
        axis(2, cex.axis=1.2, at=xat, labels=xlb)
        axis(2, at=pretty(par("usr")[3:4]), labels=FALSE)
        axis(1, at=1:length(cols), labels=FALSE)
        figlabel(paste0("n=",nrow(aam1),"  "),
                 pos="topright", font=1, cex=1.3)
        figlabel(lb, pos="topleft", cex=1.3, family=FONT, font=2)
        diffLogo_addPvals(dfop, ymin=mxo, levels=psig)
        ##mtext("encoded", 4,-.25,adj=.05)
        
        ## incorporated
        par(mai=nmai, mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
        myDiffLogo(dfnb, sparse=TRUE, ymin=0, ymax=mxn)
        xat <- pretty(par("usr")[3:4])
        xlb <- c("0",xat[length(xat)-1])
        xat <- c(0,xat[length(xat)-1])
        axis(2, cex.axis=1.2, at=xat, labels=xlb)
        axis(2, at=pretty(par("usr")[3:4]), labels=FALSE)
        ##mtext("JS divergence", 2, 1.3)
        axis(1, at=1:length(cols), labels=axlab, cex.axis=1.2)
        if ( 0 %in% cols )
            axis(1, at=which(cols==0), labels="AAS", las=1, cex.axis=1.2)
        ##text(par("usr")[2], mxn*.75, "incorporated", pos=2)
        diffLogo_addPvals(dfnp, ymin=mxn, levels=psig)
        ##mtext("incorporated", 4,-.25,adj=.05)
        dev.off()

    }
}

## y-axis label for tight motif plot
## figure heights
ht <- 3
omai <- c(.15,.25,.25,.15)
nmai <- c(.25,.25,.15,.15)
ymai <- c(.25,0,.25,0)
## heights
htn <- ht/2 + nmai[1] + nmai[3]
hto <- ht/2 + omai[1] + omai[3]
 
plotdev(file.path(mfig.path,paste0("logos_ylab")),
        height=hto+htn , width=.7, res=300, type=ftyp)
par(mai=ymai)
plot(1,1, axes=FALSE, col=NA, col.axis=NA)
text(1.25,1,label="JS divergence", xpd=TRUE, srt=90, cex=1.5)
rect(xleft=.6, ybottom=1.02, xright=1.05, ytop=par("usr")[2],
     xpd=TRUE, col="darkgray")
text(0.8,1.22,label="encoded", xpd=TRUE,srt=90, cex=1.6,
     col="white", font=2)
rect(xleft=.6, ybottom=par("usr")[3],
     xright=1.05, ytop=1.02, xpd=TRUE, col="darkgray")
text(0.82,.79,label="incorporated", xpd=TRUE,srt=90, cex=1.6, 
     col="white", font=2)
dev.off()



## RAAS dotplot motifs

tmtm <- tmtf
value <- "RAAS"

### MAP MOTIF TABLE w unique protein sites TO FULL RAAS TABLE
## THIS is used below for RAAS profiles
TIDX <- match(paste(tmtm$ensembl, tmtm$pos, tmtm$fromto),
              paste(cdat$ensembl, cdat$pos, cdat$fromto))

## is anyone missing
ina <- which(is.na(TIDX))
if ( length(ina)>0 ) 
    stop("TODO:", length(ina), "missing from unique site table.\n")

### CALCULATE ALL MOTIF RAAS PROFILES


## use RAAS profile with matrix input (new 20240604)
    
tcls <- classes[TIDX,]
ovm <- raasProfile(x=tmtm, id="SAAP", 
                   rows=tcls, cols="Dataset",
                   bg=TRUE, value=value, 
                   col.srt=uds,
                   use.test=use.test, do.plots=FALSE,
                   xlab=xl.raas,
                   verb=0)

## sort
ovm <- sortOverlaps(ovm, axis=2, srt=selected)

## SELECTED MOTIS FOR MAIN
msrt <- c(
    "KRAQ",
    "CCxCC",
    "MMxMM",
    "WWxWW"
)

ovs <- sortOverlaps(ovm, axis=2, srt=msrt)
plotProfiles(ovs,
             fname=file.path(mfig.path,paste0("motifs_",SETID,"")),
             mai=c(0.75,CMAIL,0.05,.5), ttcols=ttcols, value="median",
             p.min=p.min, p.txt=p.txt,
             dot.sze=dot.sze, p.dot=p.dot,
             rlab=LAB,  ftyp=ftyp,
             mtxt="", mtxt.line=2.3,
             vcols=acols, vbrks=abrks,
             bg=NA,  plot.all=FALSE,
             gcols=gcols, ffam=FONT)



### RAAS dotplot: CONSERVATION - DISORDER
## TODO: move this to structure.R

ovls <- raasProfile(x=tmtm, id="SAAP", 
                    rows="iupred3.bins", cols="MMSeq2.bins",
                    row.srt=c(rev(levels(iupred3.bins)),"na"),
                    col.srt=c("na", levels(MMSeq2.bins)),
                    bg=FALSE, value=value, 
                    use.test=use.test, do.plots=FALSE,
                    xlab=xl.raas,
                    verb=1)
omai <- c(.5,CMAIL,.5,.5)
nw <- ncol(ovls$p.value)*.2 + omai[2] + omai[4]
nh <- nrow(ovls$p.value)*.2 + omai[1] + omai[3]

plotdev(file.path(mfig.path,paste0("classes_conservation_disorder_raas")),
        height=nh, width=nw, res=300, type=ftyp, bg=NA)
par(mai=omai, mgp=c(1.3,.3,0), tcl=-.05, family=FONT)
dotprofile(ovls, value="median",
           p.dot=p.dot, dot.sze=dot.sze, vbrks=abrks, vcols=acols, 
           xlab=NA, ylab=NA, show.total=TRUE, tot.cex=.8)
polygon(y=c(2, nrow(ovls$p.value), nrow(ovls$p.value)),
        x=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(2, at=1, label="na", las=2)
mtext("disorder",2, 1.3)
polygon(x=c(2, ncol(ovls$p.value), ncol(ovls$p.value)),
        y=c(-.2,-.6,.2), xpd=TRUE, col="#aaaaaa", border=1)
axis(1, at=1, label="na", las=2)
mtext("conservation",1, 1.3)
dev.off()


