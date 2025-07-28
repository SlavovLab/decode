
library(readxl)
library(Biostrings) # for genetic code, blosum62, etc
library(segmenTools)
source("~/programs/segmenTools/R/coor2index.R")
source("~/programs/segmenTools/R/parsers.R")

## liftOver utility from UCSC
liftOver <- file.path("~/programs/ucsc_utils/liftOver")

## GENERATE HUMAN tRNA ANNOTATION

## TODO: load official annotation

## output
out.file <- file.path("~/data/mammary/",
                      "codons_GRCh38.tsv")

## input files
codonmap <- file.path("~/data/mammary/originalData",
                      "torres19_s3.xlsx")
codonclasses <- file.path("~/data/mammary/originalData",
                          "gingold14_mmc2.xls")
## liftOver chain
chain <-  file.path("~/data/mammary/originalData","hg19ToHg38.over.chain.gz")

## Torres et al. 2019 : alignments were performed against the human
## reference genome hg38. Human hg38 predicted tRNA genes were
## downloaded from the GtRNAdb v2.0 (January 2016) (2).

codons <- rbind(
    cbind(readxl::read_xlsx(codonmap, sheet=2), excluded=FALSE),
    cbind(readxl::read_xlsx(codonmap, sheet=1), excluded=TRUE))

## putz it
colnames(codons) <- sub("space", "chr", colnames(codons))
colnames(codons) <- tolower(colnames(codons))

codons$strand <- gsub("\\(|\\)","", codons$strand)

## assign anticodons
codons$codons <- sapply(codons$anticodon, revcomp)


## Gingold et al. 2014: tRNA classification prolif/diff/other
g14 <- as.data.frame(readxl::read_xls(codonclasses, skip=2))
g14a <- cbind(g14[1:118,1:5], class="differentiation")
g14b <- cbind(g14[1:121,7:11], class="proliferation")
g14c <- cbind(g14[1:273,13:17], class="other")
colnames(g14a) <- colnames(g14b) <- colnames(g14c) <-
    c("chr","start","end","AA","Anticodon", "class")

ccls <- rbind.data.frame(g14a,g14b,g14c)

ccls$locus <- paste0(ccls$chr,":",ccls$start, "-", ccls$end)

## ORDER BY START to allow liftOver via bed
ccls <- ccls[order(as.numeric(ccls[,"start"])), ]
ccls <- ccls[order(ccls[,"chr"]), ]

## NOTE: locations differ, tRNA set by Gingold et al.
## likely refers to an older genome version
## match(codons$locus, ccls$locus)
coor <- cbind(ccls[,c("chr", "start", "end")], strand="+")
coor$chr <- sub("chr", "", coor$chr)
coor$chr <- sub("_.*", "", coor$chr)
tmpfile <- tempfile()
## NOTE: automatically generated name corresponds to row number
## in ccls!
bed <- coor2bed(coor,  file=tmpfile)

## call liftOver and re-load coordinates
tmpfile2 <- tempfile()
tmpfile3 <- tempfile()
cmd <- paste(liftOver, tmpfile, chain, tmpfile2, tmpfile3)
system(cmd)

ncoor <- bed2coor(tmpfile2, header = c("chr", "start", "end",
                                       "name", "score","strand"))

ccls$locus.hg38 <- NA
idx <- as.numeric(sub("id","", ncoor$name))
ccls$locus.hg38[idx] <-  paste0("chr",ncoor$chr,":",ncoor$start, "-",ncoor$end)

## find locus in codon table
codons$CL.gingold14 <- ccls$class[match(codons$locus, ccls$locus.hg38)]

## fill up missing with simple anticodon match
## NOTE: ATG codon for Met is differentiated between internal and start
## so this is a bit unclear
CL.gingold14a <- ccls$class[match(codons$anticodon, ccls$Anticodon)]

missing <- which(is.na(codons$CL.gingold14))

codons$CL.gingold14[missing] <- CL.gingold14a[missing]

##codons[which(codons$CL.gingold14a!=codons$CL.gingold14),]

## MANUAL: add codons recognized by anticodon modifications

## check missing codons
GENETIC_CODE[which(!names(GENETIC_CODE)%in%codons$codons)]
## [1] "TTT" "CTC" "CAT" "CGC" "ACC" "GGT"

## The redundant genetic codons NNU and NNC (where N is A, T, G, or C)
## specify the same amino acid and are decoded by their cognate tRNAs,
## which contain either a guanosine or a modified base in the wobble
## position of the anticodons.

## tRNALeu(CAA)-specific 5-methyl-cytosine (m5C) modification in the
## wobble anticodon position enhanced translation of stress response
## genes.53

## MISSING: guesses
## Phe: TTT is decoded by GAA via T:G wobble ?? GUESS
## His: CAT is decoded by GTG via T:G wobble ?? GUESS
## Gly: GGT is decoded by GCC via T:G wobble ?? GUESS

## Leu: CTC is decoded by AAG -> IAG via I:C ?? GUESS
## Arg: CGC is decoded by ACG -> ICG via I:C ?? GUESS
## Thr: ACC is decoded by AGT -> IGT via I:C - TODO: reference

## wobble U:G
codons$codons[codons$anticodon=="GAA"] <-
    paste0(codons$codons[codons$anticodon=="GAA"],",TTT")
codons$codons[codons$anticodon=="GTG"] <-
    paste0(codons$codons[codons$anticodon=="GTG"],",CAT")
codons$codons[codons$anticodon=="GCC"] <-
    paste0(codons$codons[codons$anticodon=="GCC"],",GGT")
## wobble I:(TCA)
codons$codons[codons$anticodon=="AAG"] <-
    paste0(codons$codons[codons$anticodon=="AAG"],",CTC")
codons$codons[codons$anticodon=="ACG"] <-
    paste0(codons$codons[codons$anticodon=="ACG"],",CGC")
codons$codons[codons$anticodon=="AGT"] <-
    paste0(codons$codons[codons$anticodon=="AGT"],",ACC")

## TODO: allow all wobbles?, e.g.
## 1st anticodon position G pairs with 3rd codon position C and T,

## TODO: order and rename columns
write.table(codons, file=out.file, sep="\t", na="", quote=FALSE,
            row.names=FALSE)
