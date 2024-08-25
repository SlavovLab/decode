library(segmenTools)
options(stringsAsFactors=FALSE)

## DECODE DATA
proj.path <- file.path(Sys.getenv("DECODE")) 
if ( proj.path=="" ) # author's local path
    proj.path <- "/home/raim/data/decode"

dat.path <- file.path(proj.path,"originalData")
add.path <- file.path(proj.path,"additionalData")

fasta <- file.path(dat.path, "Homo_sapiens.GRCh38.pep.all.fa.gz")
saap.file <- file.path(add.path, "All_SAAP_protein_filter_df.txt.gz")

out.file <- file.path(proj.path,"processedData","all_proteins.fa")

if ( !file.exists(saap.file) )
    stop("input file is not available; make sure that either DECODE is set in",
         "the calling environment or that the variable proj.path in this",
         "script")

dat <- read.delim(saap.file) 



fas <- readFASTA(fasta, grepID=TRUE)
## remove version tags
names(fas) <- sub("\\.[0-9]*", "", names(fas))

## EXTRACT leading razor proteins with mutation tags

## first split by ;
## column with protein matches
PMATCH <- "Leading.razor.proteins..all.TMT." 
plst <- strsplit(dat[,PMATCH],  ";")
## clean up strings
plst <- lapply(plst, function(x) gsub("\\['","",gsub("\\']","",x)))
## split by ', '
plst <- lapply(plst, function(x) unlist(strsplit(x, "', '")))

## SPLIT BY MAPPING PIPELINE/PROTEIN ANNOTATION
unip.lst <- lapply(plst, function(x)
    unique(sub("-1","",sub("^CON__","",x[grep("^CON",x)]))))
strg.lst <- lapply(plst, function(x) x[grep("^STRG",x)])
ensp.lst <- lapply(plst, function(x) x[grep("^ENSP",x)])
    
## get unique ensembl proteins w/o mutation tags
ensm <- unique(unlist(lapply(ensp.lst,
                             function(x) grep("_",x, value=TRUE))))
    
## convert mutation tags to list
muts <- sub(".*_","", ensm)
names(muts) <- ensm # sub("_.*","", ensm)
muts <- strsplit(muts, ",")

## loop over all mutated proteins, replace mutations
## and store full sequence


nseq <- list()
errors <- matrix(FALSE, nrow=length(muts), ncol=3)
colnames(errors) <- c("nprotein","mform","merr")
for ( i in seq_along(muts) ) {
    
    sid <- sub("_.*","", names(muts)[i])
    seq <- fas[[sid]]$seq

    ft <- do.call(cbind, strsplit(muts[[i]], "[0-9]+"))
    mseq <- NULL
    if ( is.null(seq) ) {
        cat(paste("PROBLEM:", i ,"no sequence found:",sid,"\n"))
        errors[i,"nprotein"] <- TRUE
    } else if ( nrow(ft)==1 ) {
        cat(paste("PROBLEM:", i ,"wrong mutation format:",muts[[i]],
                  "\t", sid,"\n"))
        errors[i,"mform"] <- TRUE
    } else {
        ntar <- try(mutatePositions(seq, muts[[i]]))
        if ( class(ntar)=="try-error" ) {
            cat(paste("PROBLEM:", i ,"introducing mutations failed:",
                      paste(muts[[i]], collapse=";"),
                      "\t",ntar[1], "\t", sid, "\n"))
            errors[i,"merr"] <- TRUE
        } else nseq[[names(muts)[i]]] <- list(desc=names(muts)[i],
                                              details=NULL,
                                              seq=ntar)
    }  
}


cat(paste0("MUTATION REPLACEMENT SUMMARY:\n",
           "\tno sequence found: ", sum(errors[,"nprotein"]),"\n",
           "\twrong mutation format: ", sum(errors[,"mform"]),"\n",
           "\tmutation failure: ", sum(errors[,"merr"]),"\n",
           "\tsuccessful replacements: ", length(nseq),"\n"))

## fix names for blast
onames <- names(nseq)
names(nseq) <- gsub(",","_", names(nseq))

## cut names at length 47, incl. duplicate tags <50 for blastdb
nnms <- sapply(names(nseq),
               function(x)
                   paste0(unlist(strsplit(x,""))[1:min(c(47,nchar(x)))],
                          collapse=""))
nnms <- tagDuplicates(nnms)

## write out short/long name mapping
longnms <- cbind(onames, nnms)
write.table(file=sub(".fa", ".tsv", out.file), sep="\t",
            x=longnms, quote=FALSE, col.names=FALSE, row.names=FALSE)

## rename fasta
names(nseq) <- nnms

## add to fasta and write out
ofas <- append(fas,nseq)

## write out fasta
sink(file=file(out.file, open = "wt"))
for ( i in seq_along(ofas) ) 
    cat(paste0(">", names(ofas)[i], "\n", ofas[[i]]$seq, "\n"))
sink()


