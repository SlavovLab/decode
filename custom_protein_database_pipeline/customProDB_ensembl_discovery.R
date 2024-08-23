#!/usr/local/bin Rscript

args=commandArgs(trailingOnly=TRUE)
download.file("https://github.com/chambm/customProDB/archive/13cbebf9e3f26507afbbb661c9e985da2fd0840a.zip", "customProDB.zip", quiet=TRUE)
unzip("customProDB.zip")
devtools::load_all("customProDB-13cbebf9e3f26507afbbb661c9e985da2fd0840a")

library(RMariaDB)
library(rtracklayer)
library(data.table)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(Biostrings)
library(VariantAnnotation)
library(sqldf)
library(stringr)

#testdir <- '/Users/shiri/Google Drive/My Drive/Mistranslation Project/Proteogenomics_database/Testing_customProDB_062221/'
#bamFile <- paste0(testdir,'X20BR007.sorted.header.bam')
#vcfFile <- paste0(testdir,'X20BR007.vcf')
#annoPath <- paste0('/Users/shiri/Google Drive/My Drive/Mistranslation Project/genomes/ensembl/customProDB_files/')
#outputPath <- paste0(testdir,'output')

bamFile <- args[1]
vcfFile <- args[2]
annoPath <- args[3]
outputPath <- args[4]


load(paste0(annoPath, 'exon_anno.Rdata'))
load(paste0(annoPath, 'proseq.Rdata'))
load(paste0(annoPath, 'procodingseq.Rdata'))
load(paste0(annoPath, 'ids.Rdata'))
#load(paste0(annoPath, 'dbsnpinCoding.Rdata'))
txdb <- loadDb(paste0(annoPath, 'txdb.sqlite'))


calculateRPKM <- function(bamFile, exon, proteincodingonly = TRUE, ids=NULL,...){
  if(proteincodingonly&&is.null(ids)){
    stop("must supply ids mapping information if you choose proteincodingonly= TRUE")}
  anno <- GRanges(seqnames = exon$chromosome_name,ranges = IRanges(start=exon$exon_chrom_start,end=exon$exon_chrom_end), strand = exon$strand,
                  tr_name = exon$tx_name)
  
  targets <- scanBamHeader(bamFile)[[1]][['targets']]
  which <- GRanges(names(targets), IRanges(1, unname(targets)))
  all_tr <- c()
  readbychr <- c()
  
  for (i in seq_along(which)){
    param <- ScanBamParam(which=which[i], what=character())
    aln <- readGAlignments(bamFile, param=param)
    galn <- granges(aln)
    keepseqs <- seqlevels(galn)[-c(26:length(seqlevels(galn)))]
    anno <- keepSeqlevels(anno,keepseqs, pruning.mode='coarse')
    galn <- keepSeqlevels(galn,seqlevels(anno), pruning.mode='coarse')
    if(length(galn)>0){
      anno_1 <- anno[seqnames(anno)==seqnames(galn)[1]]
      exon_len <- as.data.frame(cbind(values(anno_1)[,'tr_name'],width(anno_1)))
      exonlenByTrans <- tapply(as.numeric(as.character(exon_len$V2)),exon_len$V1,sum)
      
      readcount <- countOverlaps(anno_1,galn)
      names(readcount) <-  values(anno_1)[,'tr_name']
      countbypro <- tapply(readcount, names(readcount), sum)
      RPK <- countbypro/(exonlenByTrans/1000)
      
      all_tr <- c(all_tr,RPK)    
      readbychr <- c(readbychr,sum(readcount))    
    }
  }
  
  totalReads <- sum(readbychr)
  RPKM <- all_tr/(totalReads/1e+06)

  options(stringsAsFactors=FALSE)
  
  if(proteincodingonly==TRUE){
    proex <- RPKM
    names(proex) <- ids[match(names(RPKM),ids[,'tx_name']),'pro_name']    
    proex <- proex[which((names(proex)!=''))]
    
    proex
  }else{
    RPKM
  }
}
RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
outf1 <- paste0(outputPath, 'rpkm.fasta')
print('RPKM done')

Outputproseq <- function(rpkm, cutoff="30%", proteinseq, outfile, ids, ...)
{
  if(grepl('%', cutoff)){
    cutoff <- quantile(rpkm, as.numeric(gsub('%', '', cutoff))/100)
  }else cutoff <- as.numeric(cutoff)
  s<-rpkm[rpkm >= cutoff]
  seqs <- proteinseq[proteinseq[, 'pro_name'] %in% names(s), ]
  v <- s[seqs[, 'pro_name']]
  seqs <- cbind(seqs, v)
  
  ftab <- merge(ids,seqs, by.x='pro_name', by.y='pro_name', all=FALSE, 
                stringsAsFactors=FALSE)
  ftab <- ftab[order(ftab[, 'v'], decreasing=TRUE), ]
  
  tmp <- apply(ftab, 1, function(x) 
    paste('>', x['pro_name'], " |", round(as.numeric(x['v']), 4), "|", 
          x['tx_name.x'], "|", x['gene_name'], "|", x['description'], '\n', 
          unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
  
  write(tmp, file=outfile)
  
}
Outputproseq(RPKM, 1, proteinseq, outf1, ids)
print('Proseq done')

vcf <- InputVcf(vcfFile)
idx_snv <- which(values(vcf[[1]])[['TYPE']] == "snp")
SNVvcf <- vcf[[1]][idx_snv]
idx_indel <- which((values(vcf[[1]])[['TYPE']] == "del") | (values(vcf[[1]])[['TYPE']] == "ins"))
indelvcf <- vcf[[1]][idx_indel]

#get position tables of SNVs and Indels
SNVloc <- Varlocation(SNVvcf,txdb,ids)
indelloc <- Varlocation(indelvcf,txdb,ids)
postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding=NULL, COSMIC=NULL)
postable_indel <- Positionincoding(indelvcf, exon)

#output variation table and variation protein sequence caused by SNVs
txlist <- unique(postable_snv[, 'txid'])
codingseq <- procodingseq[procodingseq$tx_id %in% txlist[[1]], ]

aaVariation <-  function(position_tab, coding, ...){
  #options(stringsAsFactors=FALSE)
  position_tab <- as.data.table(position_tab)
  old <- options(stringsAsFactors = FALSE, gsubfn.engine = "R")
  on.exit(options(old), add = TRUE)  
  
  setkey(position_tab, txid)
  coding$tx_id = as.integer(coding$tx_id)
  setDT(coding, key="tx_id")
  mtable <- unique(position_tab[coding])
  
  suppressWarnings(mtable <- sqldf("SELECT * FROM 'mtable' ORDER BY genename, txid, pincoding"))
   #mtable <- merge(position_tab, coding,by.x='txid', by.y='tx_id', all=F, 
  #                stringsAsFactors = FALSE)
   iub_mul <- list("AC"="M", "CA"="M",
                  "AG"="R", "GA"="R",
                  "AT"="W", "TA"="W",
                  "CG"="S", "GC"="S",
                  "CT"="Y", "TC"="Y",
                  "GT"="K", "TG"="K",
                  "CGT"="B", "CT,"="B", "GC,"="B", "GTC"="B", "TCG"="B", "TGC"="B",
                  "AGT"="D", "ATG"="D", "GAT"="D", "GTA"="D", "TAG"="D", "TGA"="D",
                  "ACT"="H", "ATC"="H", "CAT"="H", "CTA"="H", "TAC"="H", "TCA"="H",
                  "ACG"="V", "AGC"="V", "CAG"="V", "CGA"="V", "GAC"="V", "GCA"="V",
                  "A"="A", "T"="T", "G"="G", "C"="C")
  
  iub <- list("M"=c("A","C"),"R"=c("A","G"),"W"=c("A","T"),"S"=c("C","G"),"Y"=c("C","T"),"K"=c("G","T"),
              "B"=c("G","T","C"),"D"=c("G","T","A"),"H"=c("A","T","C"),"V"=c("G","A","C"),
              "A"="A","T"="T","G"="G","C"="C")
  
  #forward <- c("A","T","G","C","M","R","W","S","Y","K","B","D","H","V")
  #reverse <- c("T","A","C","G","K","Y","W","S","R","M","V","H","D","B")
  #nucle <- cbind(forward, reverse)
  
  index <- which(nchar(mtable$varbase) > 1)
  var_new <- unlist(lapply(mtable$varbase, function(x) if(x %in% iub_mul){iub_mul[[x]]} else{'X'}))
  
  #index <- which(nchar(as.character(mtable[, 'varbase'])) > 1)
  #var_mul <- strsplit(as.character(mtable[index, 'varbase']), ',')
  #ref_mul <- strsplit(as.character(mtable[index, 'refbase']), ',')
  
  #logi_list <- lapply(seq_along(ref_mul), function(z) lapply(strsplit(ref_mul[z][[1]], ""), function(x) lapply(strsplit(var_mul[z][[1]], ""), function(y) (x==y))))
  #var_new <- lapply(seq_along(var_mul), function(x) strsplit(var_mul[x][[1]], "")[[1]][which(logi_list[[x]][[1]][[1]]==FALSE)])
  #ref_new <- lapply(seq_along(ref_mul), function(x) strsplit(ref_mul[x][[1]], "")[[1]][which(logi_list[[x]][[1]][[1]]==FALSE)])
  
  #vars <- as.character(mtable[, 'varbase'])
  #vars[index] <- var_new
  #class(mtable[, 'varbase']) <- 'character'
  
  #refs <- as.character(mtable[, 'refbase'])
  #refs[index] <- ref_new
  #class(mtable[, 'refbase']) <- 'character'
  
  mtable$varbase <- var_new
  
  
  strand <- mtable$strand
  pincoding <- mtable$pincoding
  
  
  #strand <- as.character(mtable[,'strand'])
  #pincoding <- mtable[,'pincoding']
  
  txCodons = sqldf::sqldf(paste0("SELECT genename, txname, txid, proname, chr, strand, pos, refbase, varbase,
                           pincoding, coding,",
                                 ifelse("rsid" %in% colnames(mtable), "rsid, ", ""),
                                 ifelse("cosid" %in% colnames(mtable), "cosid, ", ""),
                                 "(ROUND((pincoding + 0.5)/3)-1)*3+1 AS CodonStart,
                           GROUP_CONCAT((CAST(pincoding AS INT) || ':' || varbase), ':') AS CodonVariants,
                           SUBSTR(coding, (ROUND((pincoding + 0.5)/3)-1)*3+1, 3) AS RefCodon
                           FROM 'mtable'
                           GROUP BY txid, (ROUND((pincoding + 0.5)/3)-1)*3+1
                           ORDER BY genename, txid, pincoding"));
  
  updateVar <- function(v, codonStart, refCodon, strand) {
    vars = unlist(stringr::str_split(v, stringr::fixed(":")))
    indices = as.numeric(vars[c(TRUE, FALSE)])-as.numeric(codonStart)+1
    varCodon = refCodon
    if (strand=="+")
      for(i in 1:length(indices))
        substr(varCodon, indices[i], indices[i]) = vars[c(FALSE, TRUE)][i]
    else
      for(i in 1:length(indices))
        substr(varCodon, indices[i], indices[i]) = fastComplement(vars[c(FALSE, TRUE)][i])
    varCodon
  }
  .complements <- c("A"="T","T"="A",
                    "G"="C","C"="G",
                    "M"="K","K"="M",
                    "R"="Y","Y"="R",
                    "W"="W","S"="S",
                    "D"="H","H"="D",
                    "B"="V","V"="B")
  
  .fastComplement = function(base) .complements[base]
  fastComplement = function(naString)
  {
    sapply(lapply(strsplit(naString, ""), .fastComplement), paste0, collapse="")
  }
  
  varcode = mapply(FUN = updateVar, txCodons$CodonVariants, txCodons$CodonStart, txCodons$RefCodon, txCodons$strand, USE.NAMES=FALSE)
  fastTranslate = function(codon) tryCatch(GENETIC_CODE[[codon]], error=function(e) "X")
  
  varaa = vector('list', length(varcode))
  vartype <- vector('character', length(varcode))
  aaref <- vector('character', length(varcode))
  aapos <- vector('integer', length(varcode))
  aavar <- vector('character', length(varcode))
  for (i in 1:length(varcode)){
    #if (show_progress) { setTxtProgressBar(pb, i) }
    aaref[[i]] = fastTranslate(txCodons$RefCodon[[i]])
    aapos[[i]] = ceiling(txCodons$CodonStart[[i]]/3)
    
    # if there are no ambiguous bases, simply do a quick translation
    if (!stringi::stri_detect_regex(varcode[[i]], "[^ACGT]")) {
      varaa[[i]] = fastTranslate(varcode[[i]])
    } else {
      # expand the ambiguous bases into combinations of unambiguous bases
      tt = iub[unlist(strsplit(varcode[[i]], split = ""))]
      combine = expand.grid(tt, KEEP.OUT.ATTRS=F, stringsAsFactors=F)
      vcodes = apply(combine, 1, paste0, collapse='')
      varaa[[i]] = paste(setdiff(unique(lapply(vcodes, fastTranslate)), aaref[[i]]), collapse=',')
    }
    
    cur_aaref = aaref[[i]]
    cur_varaa = varaa[[i]]
    
    if(is.na(match(cur_aaref, cur_varaa))) {
      vartype[[i]] = 'non-synonymous'
      aavar[[i]] = cur_varaa
    } else {
      varaaunique <- cur_varaa[-match(cur_aaref, cur_varaa)]
      if(length(varaaunique)==0) {
        vartype[[i]] = 'synonymous'
        aavar[[i]] = unique(unlist(cur_varaa))
      } else {
        vartype[[i]] = 'non-synonymous'
        aavar[[i]] = paste(varaaunique, collapse='')
      }
    }
  }
  # return input table with new columns added
  txCodons$varcode = unlist(varcode)
  txCodons$vartype = vartype
  txCodons$aaref = aaref
  txCodons$aapos = aapos
  txCodons$aavar = aavar
  as.data.frame(txCodons)
}
mtab <- aaVariation (postable_snv, codingseq)

SNV_fasta_file <- paste0(outputPath, 'SNV.fasta')
OutputVarproseq(mtab, proteinseq, SNV_fasta_file, ids)

outf_mtab <- paste0(outputPath,'SNV.tab')
write.table(mtab, file=outf_mtab, sep='\t', quote=F, row.names=F)
print('variation done')

#output fasta containing indel/frame shift proteins
txlist_indel <- unique(postable_indel[, 'txid'])
codingseq_indel <- procodingseq[procodingseq$tx_id %in% txlist_indel[[1]], ]
indel_fasta_file <- paste0(outputPath, 'indel.fasta')
Outputaberrant(postable_indel, coding=codingseq_indel, proteinseq=proteinseq, outfile=indel_fasta_file, ids=ids)
print('indel done')
