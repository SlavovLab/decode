
library(coRdon)

## read the coding sequence fasta and write a table
## of codon counts for each transcript

## codon frequency from generated coding region fasta
mam.path <- Sys.getenv("MAMDATA")

transcr <- file.path(mam.path, "processedData", "coding.fa")

out.file <- sub("\\.fa", "_codons.tsv", transcr)

## read transcripts fasta
seq <- readSet(file=transcr)
## generate codon table
cod <- codonTable(seq)

cc <- codonCounts(cod)
rownames(cc) <- sub(",.*", "", names(seq))

## write out global codon frequencies
out <- cbind.data.frame(ID=rownames(cc), cc)
write.table(out, file=out.file, sep="\t", row.names=FALSE, quote=FALSE)


