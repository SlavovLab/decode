library(stringr)
library(dplyr)
library(ggplot2)
library(seqinr)
library(protViz)
library(MSnbase)

Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path,set.attributes = T,whole.header = T)
  convert_mouse <- names(convert_mouse)
  parse_row<-grep("GN=",convert_mouse, fixed=T)
  split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
  gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
  prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
  prot_parse <- grep("|",prot, fixed=T)
  gene_parse <- grep(" ",gene, fixed=T)
  split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
  split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
  split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
  split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))
  
  return(convert_mouse)
}

# Download supplemental table 2
SAAP <- read.csv('/Users/andrewleduc/Downloads/Supplemental_Data_2.SAAP_proteins.csv')



# Set path to Analysis Dependencies/Pull_down_analysis/ folder 
path <- '/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/SuppData/2024_Tsour/Pipline_output/analysis_dependencies/Pull_down_analysis/'
# Searched pull down data from 
raw <- read.delim(paste0(path,'PullDown_search.tsv'),sep = '\t')

raw <- raw %>% filter(Intensity != 0)
SAAP <- SAAP %>% filter(is.na(Mean.precursor.RAAS)==F)

raw$Spectrum.File <- str_remove(raw$Spectrum.File,'.pep.xml')
raw$Spectrum.File <- str_remove(raw$Spectrum.File,'_rank1')
raw$Spectrum.File <- str_remove(raw$Spectrum.File,'_rank2')
raw$Spectrum.File <- str_remove(raw$Spectrum.File,'_rank3')
raw$Spectrum.File <- str_remove(raw$Spectrum.File,'_rank4')
raw$Spectrum.File <- str_remove(raw$Spectrum.File,'_rank5')
raw$Spectrum.File <- str_extract(raw$Spectrum.File, "[^_]+$")

unique(SAAP_check$Uniprot_ID)
SAAP_check <- SAAP %>% filter(Uniprot_ID %in% convert_2$Uniprot)

SAAP <- SAAP %>% filter(BP %in% raw$Peptide)
SAAP <- SAAP %>% filter(SAAP %in% raw$Peptide)


df <- matrix(data = NA,nrow = 5000,ncol = 9)
colnames(df) = c('BP','SAAP','Uniprot','Pull_down_gene','RAAS','RAAS_Shiri','numb_raw_files','SAAP_PEP','BP_PEP')
df <- as.data.frame(df)

count = 0
for(i in 1:nrow(SAAP)){
  raw_hold <- raw %>% filter(Peptide %in% c(SAAP$SAAP[i],SAAP$BP[i]))

  for(j in unique(raw_hold$Spectrum.File)){
    raw_hold_hold <- raw_hold %>% filter(Spectrum.File == j)
    
    for(k in unique(raw_hold_hold$Spectrum.File)){
      
      raw_hold_hold_hold <- raw_hold_hold %>% filter(Spectrum.File == k)
      
      if(length(unique(raw_hold_hold_hold$Peptide))==2){
        count = count +1
        
        
        bp_df <- raw_hold_hold_hold %>% filter(Peptide == SAAP$BP[i])
        dp_df <- raw_hold_hold_hold %>% filter(Peptide == SAAP$SAAP[i])
        
        df$BP[count] <- SAAP$BP[i]
        df$SAAP[count] <- SAAP$SAAP[i]
        df$Uniprot[count] <- bp_df$Protein.ID[1]
        df$Pull_down_gene[count] <- raw_hold_hold_hold$Spectrum.File[1]
        
        df$RAAS[count] <- log10(max(dp_df$Intensity,na.rm = T)/max(bp_df$Intensity,na.rm = T))
        
        df$RAAS_Shiri[count] <- SAAP$Mean.precursor.RAAS[i]
        #df$numb_raw_files
        
        df$SAAP_PEP[count] <- min(1-dp_df$Probability,na.rm=T)
        df$BP_PEP[count] <- min(1-bp_df$Probability,na.rm=T)
        
      }
    }
  }
}

df <- df %>% filter(is.na(BP) == F)


plot(df$RAAS,df$RAAS_Shiri)
cor(df$RAAS,df$RAAS_Shiri)

cor.test(df$RAAS,df$RAAS_Shiri)



df_collapsed <- df %>%
  dplyr::group_by(BP, SAAP, Pull_down_gene) %>%
  dplyr::summarize(
    Uniprot        = first(Uniprot),
    RAAS           = mean(RAAS),
    RAAS_Shiri     = first(RAAS_Shiri),
    SAAP_PEP       = first(SAAP_PEP),
    BP_PEP         = first(BP_PEP),
    .groups = "drop"
  )


plot(df_collapsed$RAAS,df_collapsed$RAAS_Shiri,main = 'Pearson = 0.35')
cor(df_collapsed$RAAS,df_collapsed$RAAS_Shiri,method = 'spearman')
cor.test(df_collapsed$RAAS,df_collapsed$RAAS_Shiri)

# Human fasta download for gene name mapping from Analysis Dependancies
convert = Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Human.fasta')
convert = convert %>% filter(split_prot %in% df_collapsed$Uniprot)
colnames(convert)[1] <- 'Uniprot'
colnames(convert)[2] <- 'Gene'

df_collapsed <- df_collapsed %>% left_join(convert, by = c('Uniprot'))

convert_2 <- convert %>% filter(Gene %in% raw$Spectrum.File)

df_collapsed3 <- df_collapsed %>% filter(SAAP_PEP < .01)
df_collapsed3 <- df_collapsed3 %>% filter(BP_PEP < .01)


df_collapsed3$Matched <- df_collapsed3$Gene == df_collapsed3$Pull_down_gene
df_collapsed3$Matched[is.na(df_collapsed3$Matched)] <- F
ggplot(df_collapsed3,aes(x = RAAS,y = RAAS_Shiri,color = Matched))+ geom_point(size=2)+
  theme_classic(base_size = 18)+ geom_abline(intercept = 0,slope = 1)+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_blank()  
  ) + scale_color_manual(values = c('grey30','darkorange')) + xlim(c(-4,2))+ylim(c(-4,2))+
  xlab('RAAS Pull down') + ylab('RAAS CPTAC')



write.csv(df_collapsed3,'~/Desktop/Supplementary_data_table_x.csv')



#### Plotting raw spectrum



# look for correct scan numbers to read data from mzml
HSPD1 <- read.delim('/HSPD1.tsv')
HSPD1_BP <- HSPD1 %>% filter(peptide == 'IMQSSSEVGYDAMAGMFVNMVEK') # Scan number 18951
HSPD1_DP <- HSPD1 %>% filter(peptide == 'IMQSSSEVGYDAMAGDFVNMVEK') # Scan number 16257


prosit <- read.csv('prosit_koinapy_predictions_pulldowns.csv')

## Converted MZML
ms <- readMSData("HSPD1.mzML", mode = "onDisk")


##############################################################################################################
# Plot spectra for DP


# Suppose your peptide sequence is:
pep <- "IMQSSSEVGYDAMAGDFVNMVEK"   # replace with your own sequence

# Calculate all b‐ and y‐ion masses (singly charged by default):
frags <- fragmentIon(pep)

frags <- frags[[1]]



spec <- ms[ acquisitionNum(ms) == 16257 ][[1]]   # scan 16258  (ms level 2)
#plot(spec)  



df_pep_plot_BP <- data.frame(theo = c(round(frags$b,1), round(frags$y,1)), anno = c(paste0('b',1:23), paste0('y',1:23)),type = c(rep('b',23),rep('y',23)))   

df_pep_plot_BP <- df_pep_plot_BP %>% filter(theo %in% round(spec@mz,1))
df_real <- data.frame(mz = round(spec@mz,1),intensity = spec@intensity)

df_real <- df_real %>% filter(mz %in% df_pep_plot_BP$theo)


df_real <- df_real %>% distinct(mz,.keep_all = T)

df_real <- df_real[order(df_real$mz),]
df_pep_plot_BP <- df_pep_plot_BP[order(df_pep_plot_BP$theo),]

df_pep_plot_BP$Intensity = df_real$intensity

df_pep_plot_BP$Intensity <- df_pep_plot_BP$Intensity/max(df_pep_plot_BP$Intensity)


df_pep_plot_BP$anno <- NA
df_pep_plot_BP$type <- NULL


prosit_one <- prosit %>% filter(peptide_sequences == 'IMQSSSEVGYDAMAGDFVNMVEK')

prosit_one <- prosit_one %>% dplyr::select(mz,annotation,intensities)
colnames(prosit_one) <- colnames(df_pep_plot_BP)
prosit_one$Intensity <- -prosit_one$Intensity

df_pep_plot_BP_p <- rbind(df_pep_plot_BP,prosit_one)
                  


df_pep_plot_BP_p$anno <- str_sub(df_pep_plot_BP_p$anno,3,-1)
df_pep_plot_BP_p$anno <- str_sub(df_pep_plot_BP_p$anno,1,-4)


df_pep_plot_BP_p$label_peak <- T
## --- 2.  Mirror plot with annotations  --------------------------------

df_pep_plot_BP_p$spect <- 'Empirical'
df_pep_plot_BP_p$spect[df_pep_plot_BP_p$Intensity >0] <- 'Prosit'

ggplot(df_pep_plot_BP_p, aes(x = theo, xend = theo,
                             y = 0,   yend = Intensity,color = spect)) +
  geom_segment(linewidth = 0.6) +
  ## annotation layer
  geom_text(
    data = subset(df_pep_plot_BP_p, label_peak),
    aes(x = theo,
        y = Intensity + 0.05 * sign(Intensity) * max(abs(Intensity)),  # nudges label
        label = anno),
    angle = 90,
    hjust = ifelse(df_pep_plot_BP_p$Intensity > 0, 0, 1),             # above vs below
    size  = 3
  ) +
  geom_hline(yintercept = 0) +
  labs(x = "m/z", y = "Relative intensity") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_blank()  
  ) + scale_color_manual(values = c('red','black')) + xlim(c(100,1500))





##############################################################################################################
# Plot spectra for BP

# Suppose your peptide sequence is:
pep <- "IMQSSSEVGYDAMAGMFVNMVEK"   # replace with your own sequence

# Calculate all b‐ and y‐ion masses (singly charged by default):
frags <- fragmentIon(pep)

frags <- frags[[1]]



spec <- ms[ acquisitionNum(ms) == 18951 ][[1]]   # scan 16258  (ms level 2)
#plot(spec)  



df_pep_plot_BP <- data.frame(theo = c(round(frags$b,1), round(frags$y,1)), anno = c(paste0('b',1:23), paste0('y',1:23)),type = c(rep('b',23),rep('y',23)))   

df_pep_plot_BP <- df_pep_plot_BP %>% filter(theo %in% round(spec@mz,1))
df_real <- data.frame(mz = round(spec@mz,1),intensity = spec@intensity)

df_real <- df_real %>% filter(mz %in% df_pep_plot_BP$theo)


df_real <- df_real %>% distinct(mz,.keep_all = T)

df_real <- df_real[order(df_real$mz),]
df_pep_plot_BP <- df_pep_plot_BP[order(df_pep_plot_BP$theo),]

df_pep_plot_BP$Intensity = df_real$intensity

df_pep_plot_BP$Intensity <- df_pep_plot_BP$Intensity/max(df_pep_plot_BP$Intensity)


df_pep_plot_BP$anno <- NA
df_pep_plot_BP$type <- NULL


prosit_one <- prosit %>% filter(peptide_sequences == 'IMQSSSEVGYDAMAGMFVNMVEK')

prosit_one <- prosit_one %>% dplyr::select(mz,annotation,intensities)
colnames(prosit_one) <- colnames(df_pep_plot_BP)
prosit_one$Intensity <- -prosit_one$Intensity

df_pep_plot_BP_p <- rbind(df_pep_plot_BP,prosit_one)



df_pep_plot_BP_p$anno <- str_sub(df_pep_plot_BP_p$anno,3,-1)
df_pep_plot_BP_p$anno <- str_sub(df_pep_plot_BP_p$anno,1,-4)


df_pep_plot_BP_p$label_peak <- T
## --- 2.  Mirror plot with annotations  --------------------------------

df_pep_plot_BP_p$spect <- 'Empirical'
df_pep_plot_BP_p$spect[df_pep_plot_BP_p$Intensity >0] <- 'Prosit'

ggplot(df_pep_plot_BP_p, aes(x = theo, xend = theo,
                             y = 0,   yend = Intensity,color = spect)) +
  geom_segment(linewidth = 0.6) +
  ## annotation layer
  geom_text(
    data = subset(df_pep_plot_BP_p, label_peak),
    aes(x = theo,
        y = Intensity + 0.05 * sign(Intensity) * max(abs(Intensity)),  # nudges label
        label = anno),
    angle = 90,
    hjust = ifelse(df_pep_plot_BP_p$Intensity > 0, 0, 1),             # above vs below
    size  = 3
  ) +
  geom_hline(yintercept = 0) +
  labs(x = "m/z", y = "Relative intensity") +
  theme_classic()+
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_blank()  
  ) + scale_color_manual(values = c('red','black')) + xlim(c(100,1500))


















