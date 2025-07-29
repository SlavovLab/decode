# Read in text file containing all the relevant raw data downloaded from DESCRIBEPROT at the link below

## https://biomine.cs.vcu.edu/servers/DESCRIBEPROT/download_database/9606_database.json

# Path to file
mam.path <- Sys.getenv("MAMDATA")
in.file <-  file.path(mam.path,  'originalData', '9606_database.json')


out.path <- file.path(mam.path, 'processedData', 'describePROT')
dir.create(out.path)

# Text mine file line by line and extract saved scores for each amino acid
con <- file(in.file, "r")
all_line <- readLines(con)
close(con)

prot_list <- c()
conserve_scores <- c()
surface_acc <- c()
prot_diso <- c()
prot_binding <- c()
PTM <- c()
seq_p <- c()
disorder_res <- c()

flag = T
for(i in 1:length(all_line)){
  if(length(grep('"ACC\"',all_line[i])==1)==1){
    prot <- gsub('.*"(\\w+)"', '\\1', all_line[i])
    prot <- substr(prot, 1, nchar(prot) - 2)
    if(flag == T){
      prot_list <- c(prot_list,prot)
      flag = F
    }
    
  }

  # Searching each line to see if it contains information from one of the desired scores
  
  if(length(grep('MMseq2_conservation_score',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 33, nchar(all_line[i]) - 3)
    conserve_scores <- c(conserve_scores,score_save)
    flag = T
  }
  
  if(length(grep('ASAquick_normscore',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 26, nchar(all_line[i]) - 3)
    surface_acc <- c(surface_acc,score_save)
    flag = T
  }
  
  if(length(grep('DisoPROscore',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 20, nchar(all_line[i]) - 3)
    prot_diso <- c(prot_diso,score_save)
    flag = T
  }
  
  if(length(grep('SCRIBERscore',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 20, nchar(all_line[i]) - 3)
    prot_binding <- c(prot_binding,score_save)
    flag = T
  }
  
  if(length(grep('PTMbinary',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 17, nchar(all_line[i]) - 3)
    PTM <- c(PTM,score_save)
    flag = T
  }
  
  if(length(grep('flDPnn_score',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 20, nchar(all_line[i]) - 3)
    disorder_res <- c(disorder_res,score_save)
    flag = T
  }
  
  if(length(grep('"seq":',all_line[i])==1)==1){
    score_save <- substr(all_line[i], 11, nchar(all_line[i]) - 3)
    seq_p <- c(seq_p,score_save)
    flag = T
  }
  
  
}

df <- as.data.frame(prot_list)
df$conserve <- conserve_scores
df$surf_acc <- surface_acc
df$protDiso <- prot_diso
df$ProtBind <- prot_binding
df$diso_res <- disorder_res
df$seq <- seq_p


# Unlist the scores and save a seperate file for each protein which has all the scores by Amino acid residue


for(i in 1:nrow(df)){
  
  pep <- unlist(str_split(df$seq[i], ''))
  
  df_temp <- as.data.frame(pep)
  if(length(unlist(str_split(df$conserve[i], ',')))== nrow(df_temp)){
    df_temp$conserve <- unlist(str_split(df$conserve[i], ','))
  }else{
    df_temp$conserve <- rep(NA,nrow(df_temp))
  }
  
  if(length(unlist(str_split(df$surf_acc[i], ',')))== nrow(df_temp)){
    df_temp$surf_acc <- unlist(str_split(df$surf_acc[i], ','))
  }else{
    df_temp$surf_acc <- rep(NA,nrow(df_temp))
  }
  
  if(length(unlist(str_split(df$protDiso[i], ',')))== nrow(df_temp)){
    df_temp$protDiso <- unlist(str_split(df$protDiso[i], ','))
  }else{
    df_temp$protDiso <- rep(NA,nrow(df_temp))
  }
  
  if(length(unlist(str_split(df$ProtBind[i], ',')))== nrow(df_temp)){
    df_temp$ProtBind <-unlist(str_split(df$ProtBind[i], ','))
  }else{
    df_temp$ProtBind <- rep(NA,nrow(df_temp))
  }
  
  
  
  if(length(unlist(str_split(df$diso_res[i], ',')))== nrow(df_temp)){
    df_temp$diso_res <-unlist(str_split(df$diso_res[i], ','))
  }else{
    df_temp$diso_res <- rep(NA,nrow(df_temp))
  }
  
  
    write.table(df_temp, file.path(out.path, paste0(df$prot_list[i],".tsv")),
                sep = "\t", row.names = FALSE, col.names = FALSE)
  
}
