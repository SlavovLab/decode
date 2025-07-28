
library(coRdon)

## read the coding sequence fasta and write a table
## of codon counts for each transcript

## codon frequency from generated coding region fasta
mam.path <- Sys.getenv("MAMDATA")
transcr <- file.path(mam.path,"processedData","coding.fa")#"Homo_sapiens.GRCh38.
out.file <- sub("\\.fa", "_codons.tsv", transcr)
seq <- readSet(file=transcr)
cod <- codonTable(seq)

cc <- codonCounts(cod)
rownames(cc) <- sub(",.*", "", names(seq))

out <- cbind.data.frame(ID=rownames(cc), cc)
write.table(out, file=out.file, sep="\t", row.names=FALSE, quote=FALSE)

## AA frequency


if ( FALSE ) {

    library(segmenTools)
    ## test coRdon: 

    dnaHD59 <- readSet(
        file="https://raw.githubusercontent.com/BioinfoHR/coRdon-examples/master/HD59.fasta"
    )
    HD59 <- codonTable(dnaHD59)
    
    xlab <- "MILC distance from sample centroid"
    ylab <- "MILC distance from ribosomal genes"
    milc <- MILC(HD59, ribosomal = TRUE)
    
    dense2d(milc[,1], milc[,2])

    genes <- getKO(HD59)[getlen(HD59) > 80]
    
    subset <- list(half=c(rep(TRUE,20), rep(FALSE,length(seq)-20)))
    
    milc <- MILC(cod, self = FALSE,
                 subsets = subset)
}
