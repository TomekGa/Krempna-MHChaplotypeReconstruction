### HAPLOTYPE RECONSTRUCTION ###
### Tomasz Gaczorek ###
### tomasz.gaczorek@doctoral.uj.edu.pl ###
### 05.01.2023 ###

### INTRO ####
#STEP 1 - reading raw data, plotting and resolving maternal haplotypes
#REQUIRE: binary_dt.csv 
#SOFTWARE: build_tree() function requires locally installed mafft (https://mafft.cbrc.jp/alignment/software/)

### IMPORTANT ###
#The code below was built on the Linux platform. You may encounter some problems while using it on Windows.
#Known problems appear in the "TREE FOR EXON" section utilizing align_mafft() and build_tree() functions:
# align_mafft() - it calls locally installed mafft software and uses command line to align sequences. If needed, install Windows version or use online mafft version to align your sequences.
# build_tree() - it uses phangorn::bootstrap.pml() which can use multiple cores on the Linux platform. Such option is not implemented on Windows. If needed, change the argument "multi" to FALSE.

#setwd("~/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG") #setting working directory
setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG") #windows' version
library(tidyverse)
library(ggtext)
#
### INITIAL DATA PREPARATION  - NOT NEEDED IN THE FINAL VERSION #### 
# pgr <- read.csv("raw_data/pedigree.csv")
# mothers_ids <- pgr %>% filter(is.na(ID_MOTHER)) %>% pull(ID) #rows with NA means mothers
# markers <- sort(c("Iex2", "Iex3", "IIex2"))
# suff <- c("Aug_22", "Nov_22")
# 
# allele_seqs <- rbind_forMe(paths = paste0("./raw_data/",markers, "_alleles.txt"),fun_var = read_tsv) %>%
#   setNames(c("ALLELE","SEQUENCE"))
# write.csv(allele_seqs,"./generated_data/allele_seqs.csv",row.names = F,fileEncoding = "UTF-8")
# 
# binary_dt <- rbind_forMe(paths = as.vector(sapply(X = paste0("./raw_data/", markers, "_binary_genotypes_"),FUN = paste0,suff,".txt")),
#                          fun_var = read_longformat_forBinary,second_arg = sort(rep(markers,2))) %>%
#   mutate("PEDIGREE" = case_when(ID %in% mothers_ids ~ "mother",T ~ "offspring")) %>%
#   left_join(pgr,by = "ID") %>%
#   mutate("ID" = as.character(ID)) %>%
#   left_join(is_in_mother_faster(.),by = c("ID","ALLELE")) %>%
#   left_join(is_in_offspring(.),by = c("ID","ALLELE")) %>%
#   left_join(allele_seqs,by = "ALLELE") 
# write.csv(binary_dt,"generated_data/binary_dt.csv",row.names = F,fileEncoding = "UTF-8")
# #
### DATA IMPORT ####
binary_dt <- read.csv("generated_data/binary_dt.csv",stringsAsFactors = F)
filtered_binary_dt <- read.csv("generated_data/binary_dt_filtered.csv") # file generated in "FILTERING AFTER VISUAL INVESTIGATION" section below
#
### FUNCTIONS ####
# rbind_forMe <- function(paths,fun_var,second_arg = NULL){ #transform your files using provided data and bind them by rows 
#   out_list <- list()
#   if(is.null(second_arg)==T){
#     for(i in seq_along(paths)){
#       out_list[[i]] <- fun_var(paths[i])
#     }
#   } else {
#     for(i in seq_along(paths)){
#       out_list[[i]] <- fun_var(paths[i],second_arg[i])
#     }
#   }
#   as.data.frame(do.call(rbind,out_list))
# }

# read_longformat_forBinary <- function(path,marker){
#   read_tsv(path) %>% pivot_longer(2:ncol(.), names_to = "ALLELE", values_to = "PRESENT") %>% 
#     filter(PRESENT == 1) %>% select(-PRESENT) %>% mutate("MARKER" = marker)
# }

is_in_mother_faster <- function(dt){
  out <- list()
  for(i in unique(dt$ID_MOTHER)){
    if(is.na(i) == F){
      children <- dt %>% filter(ID_MOTHER == i)
      mothers_alleles <- dt %>% filter(ID == i) %>% pull(ALLELE)
      for(j in 1:nrow(children)){
        out[[length(out)+1]] <- c(children$ID[j],children$ALLELE[j],children$ALLELE[j] %in% mothers_alleles)
      }
    }
  }
  setNames(as.data.frame(do.call(rbind,out)),nm = c("ID","ALLELE","IS_IN_MOTHER"))
}

is_in_offspring <- function(dt){
  out <- list()
  for(i in unique(dt$ID_MOTHER)){
    if(is.na(i) == F){
      children_alleles <- dt %>% filter(ID_MOTHER == i) %>% pull(ALLELE) %>% unique()
      mother <- dt %>% filter(ID == i)
      for(j in 1:nrow(mother)){
        out[[length(out)+1]] <- c(mother$ID[j],mother$ALLELE[j],mother$ALLELE[j] %in% children_alleles)
      }
    }
  }
  setNames(as.data.frame(do.call(rbind,out)),nm = c("ID","ALLELE","IS_IN_OFFSPRING"))
}

clustered_order_m_or_f <- function(dt,parent = "M"){
  if(parent == "M"){
  x <- dt %>% mutate("GENOTYPE" = 1) %>% filter(IS_IN_MOTHER == T) %>%
    select(ID,ALLELE,GENOTYPE) %>%
    pivot_wider(names_from = ALLELE,values_from = GENOTYPE,values_fill = 0) #restrain to maternal alleles and change data into wider format
  } else if(parent == "F"){
    x <- dt %>% mutate("GENOTYPE" = 1) %>% ungroup() %>% 
      add_row("ID" = unique(.$ID),"IS_IN_MOTHER" = F,"ALLELE" = "alleleX","GENOTYPE" = 1) %>% #dealing with individuals lacking father's alleles
      filter(IS_IN_MOTHER == F) %>% #restrain to paternal alleles
      select(ID,ALLELE,GENOTYPE) %>%
      pivot_wider(names_from = ALLELE,values_from = GENOTYPE,values_fill = 0) #change data into wider format
  }
  d <- as.matrix(dist(x[,-1],method = "binary")) #calculate binary distance between individuals
  diag(d) <- NA
  rownames(d) <- x$ID; colnames(d) <- x$ID
  sampled_ind <- as.character(sample(x$ID,1)) #start from random individual
  IDs <- as.character(x$ID)
  for(i in 1:(ncol(d)-1)){
    minimum <- which.min(d[,sampled_ind[i]]) #index of minimum distance
    sampled_ind[i+1] <- IDs[minimum] #save individual with corresponding index
    IDs <- IDs[IDs != sampled_ind[i]] #delete individual chosen above
    d <- d[rownames(d) != sampled_ind[i],colnames(d) != sampled_ind[i]] #delete individual chosen above from dist matrix
  }
  sampled_ind
} # returns a vector of IDs sorted depending on maternal or paternal alleles similarity

align_mafft <- function(in_file,out_file,rev_comp = T,path_mafft = "/home/tomek/miniconda3/bin/mafft",method_v = "localpair"){
  initial_bin <- ape::read.FASTA(in_file)
  if(rev_comp == T){
    aligned <- ips::mafft(initial_bin,exec = path_mafft,op = 5,ep = 3,options = c("--adjustdirection"),method = method_v)
  } else {
    aligned <- ips::mafft(initial_bin,exec = path_mafft,op = 5,ep = 3,method = method_v)
  }
  ape::write.FASTA(aligned,out_file)
  aligned
} #alligning using local mafft on LINUX

build_tree <- function(alligned_DNAbin,tree_name,multi = T,boot = T){
  require(phangorn)
  phyD <- phyDat(alligned_DNAbin)
  dm <- dist.ml(phyD,model = "JC69",exclude = "pairwise")
  nj <- NJ(dm)
  fit <- pml(nj,phyD)
  if(boot == T){
    fitJC <- optim.pml(fit, model = "JC", rearrangement = "ratchet")
    bs <- bootstrap.pml(fitJC, bs=500, optNni=TRUE, multicore=multi,mc.cores = 10,control = pml.control(trace = 0))
    consensus_tree <- plotBS(midpoint(fitJC$tree), bs, type="none",p = -1)
    ape::write.tree(consensus_tree,tree_name)
    consensus_tree
  } else {
    fitJC <- optim.pml(fit, model = "JC", rearrangement = "none")
    ape::write.tree(fitJC$tree,tree_name)
    fitJC$tree
  }
}

get_union_haplotype <- function(my_list){ #intersection between individuals within a cluster + tolerance of a single allele if absent in only one individual 
  united <- Reduce(intersect,my_list)
  if(length(my_list) > 2){
    others <- unlist(my_list) %>% .[!(. %in% united)] %>% unique()
    absent <- c()
    for(i in others){
      counter <- 0
      for(j in my_list){
        if(i %in% j == F){
          counter <- counter+1
        }
      }
      absent <- c(absent,counter)
    }
    if(sum(absent == 1) > 0){
      united <- c(united,others[absent == 1])
    }
  }
  united
}
#
### PLOTTING PER MOTHER - PRIOR TO FILTERING ####
dir.create("OUTCOMES")
marker_var <- c("Iex2","Iex3")

# mother_haps <- read.csv("Mothers_haplotypes_intersect.csv") %>% group_by(ID_MOTHER,ALLELE) %>%
#   mutate("OCCURRENCE" = n()) %>% ungroup() %>%
#   mutate("HAPLOTYPE" = case_when(OCCURRENCE == 2 ~ "both",T ~ as.character(HAPLOTYPE))) %>%
#   select(-OCCURRENCE) %>% distinct()

dt_org <- binary_dt %>% 
  #filter(!(ID %in% filtering_inds)) %>% #filtering out alleles manually detected by Wiesiek
  #filter(!(ALLELE %in% filtering_alleles)) %>% #filtering out alleles manually detected by Wiesiek
  #mutate("ID" = factor(as.character(ID),levels = sort(unique(.$ID),decreasing = T))) %>% #sorting by ID
  filter(PEDIGREE == "offspring" & MARKER %in% marker_var) %>%
  select(ID,ALLELE,ID_MOTHER,IS_IN_MOTHER,MARKER) %>%
  group_by(ALLELE) %>% mutate("N_IND" = n_distinct(ID)) #%>% #number of all individuals possessing a given allele
  #left_join(mother_haps,by = c("ID_MOTHER","ALLELE")) %>% #joining mothers haplotypes
  #mutate("HAPLOTYPE" = case_when(is.na(HAPLOTYPE) ~ "none",T ~ HAPLOTYPE)) #give label to paternal alleles

for(i in unique(dt_org$ID_MOTHER)){ #for all mothers
  p_list <- list() #list for plots
  dt_m <- dt_org %>% filter(ID_MOTHER == i)
 
   for(j in c("M","F")){ #grouped by mother or father
    dt <- dt_m %>% mutate("ID" = factor(as.character(ID),levels = clustered_order_m_or_f(.,parent = j))) #sorts individuals depending on maternal or paternal alleles similarity
     
    dt <- dt %>% arrange(desc(IS_IN_MOTHER),MARKER,desc(N_IND)) %>% #sorting alleles: maternal first then most frequent
      mutate("ALLELE" = paste0(ALLELE," (",N_IND,")")) %>% #adding the number of possessors to allele labels
      # mutate("ALLELE" = case_when(HAPLOTYPE == "both" ~ paste0("<i style='color:blue4'>",ALLELE,"</i>"),
      #                             HAPLOTYPE == "1" ~ paste0("<i style='color:green'>",ALLELE,"</i>"),
      #                             HAPLOTYPE == "2" ~ paste0("<i style='color:red'>",ALLELE,"</i>"),
      #                             HAPLOTYPE == "none" ~ paste0("<i style='color:black'>",ALLELE,"</i>"))) %>% #adding color to x-labels
      mutate("ALLELE" = factor(ALLELE,levels = unique(.$ALLELE))) #changing allele variable back to factor
    
    p_list[[length(p_list)+1]] <- ggplot(dt)+
      geom_tile(aes(x = ALLELE,y = ID,fill = IS_IN_MOTHER,alpha = MARKER),col = "black")+
      theme_bw()+
      scale_fill_manual(values = c("TRUE" = "darkgreen","FALSE" = "red4"))+
      scale_alpha_manual(values = c("Iex2" = 1,"Iex3" = 0.6))+
      theme(axis.text.x = element_markdown(angle = 45,hjust = 1),panel.grid = element_blank())+
      labs(title = paste("MOTHER",i,"-",paste(marker_var,collapse = "_")), y = "OFFSPRING ID", 
           x = "ALLELE (POSSESSORS IN WHOLE DATASET)",fill = "IS IN MOTHER?")
  }
  p_final <- ggpubr::ggarrange(p_list[[1]],p_list[[2]],ncol = 1)
  
  ggsave(filename = paste0("./OUTCOMES/MY_PLOTS/prior_to_filtering/mother_",i,"_",paste(marker_var,collapse = "_"),".png")
         ,device = "png",plot = p_final,width = 20,height = 14)
}
#
### FILTERING AFTER VISUAL INVESTIGATION ####
filtering_alleles <- read.csv("generated_data/unstable.csv") %>% pull(allele)
filtering_inds <- c("24618", "24479", "24726")
filtered_binary_dt <- binary_dt %>%
  filter(!(ID %in% filtering_inds)) %>% #filtering out manually detected alleles
  filter(!(ALLELE %in% filtering_alleles)) %>% #filtering out manually detected alleles
write.csv(filtered_binary_dt,"generated_data/binary_dt_filtered.csv",row.names = F,fileEncoding = "UTF-8")
#  
### TREE FOR EXON ####
dir.create("OUTCOMES/fastas")
dir.create("OUTCOMES/trees")
#setting files' names
marker_var <- "Iex2"
fasta_name <- paste0("./OUTCOMES/fastas/",marker_var,".fasta")
fasta_aligned_name <- paste0("./OUTCOMES/fastas/",marker_var,"_aligned.fasta")
tree_name <- paste0("./OUTCOMES/trees/",marker_var,"_tree.new")
#datasets
dt <- binary_dt %>% filter(PEDIGREE == "offspring" & MARKER == marker_var)
#writting fasta
seqs <- select(dt,ALLELE,SEQUENCE) %>% distinct() %>% arrange(ALLELE)
seqinr::write.fasta(as.list(seqs$SEQUENCE),names = seqs$ALLELE,file.out = fasta_name,as.string = T)
#meta info containing info gathered during visual investigation
meta <- left_join(seqs,filtered_binary_dt,by = "ALLELE")
#alligning fasta
align_mafft(in_file = fasta_name,out_file = fasta_aligned_name)
#building tree
alligned <- ape::read.FASTA(fasta_aligned_name)
tree <- build_tree(alligned,tree_name,multi = T,boot = T)
tree$node.label[is.na(tree$node.label)] <- 0 
# plotting
library(ggtree)
plot <- ggtree(tree) %<+% meta +
  geom_tiplab(hjust = -0.2,size = 2)+ 
  geom_nodelab(geom = "label",hjust = 0.9,size = 1.5) + 
  geom_tippoint(aes(x =x + 0.01*max(x),color = unst_fr_fam),size = 1.5) + 
  scale_color_gradient(low = "green",high = "red",na.value = "transparent") + 
  theme(legend.key.size = unit(1,"cm"),legend.text = element_text(size = 5))+
  xlim(0,0.45)+
  geom_treescale(y = -1,linesize = 1,fontsize = 1.5,width = 0.05)
ggsave(paste0("./trees/",marker_var,"_tree.png"),plot,device = "png",scale = 1,limitsize = F,height = 18) #choose optimal graphical parameters
#
### PLOTTING PER MOTHER - FILTERED DATASET ####
marker_var <- c("Iex2","Iex3")

dt_org <-filtered_binary_dt %>% 
  filter(PEDIGREE == "offspring" & MARKER %in% marker_var) %>%
  select(ID,ALLELE,ID_MOTHER,IS_IN_MOTHER,MARKER) %>%
  group_by(ALLELE) %>% mutate("N_IND" = n_distinct(ID)) #number of all individuals possessing a given allele

for(i in unique(dt_org$ID_MOTHER)){ #for all mothers
  p_list <- list() #list for plots
  dt_m <- dt_org %>% filter(ID_MOTHER == i)
  
  for(j in c("M","F")){ #grouped by mother or father
    dt <- dt_m %>% mutate("ID" = factor(as.character(ID),levels = clustered_order_m_or_f(.,parent = j))) #sorts individuals depending on maternal or paternal alleles similarity
    
    dt <- dt %>% arrange(desc(IS_IN_MOTHER),MARKER,desc(N_IND)) %>% #sorting alleles: maternal first then most frequent
      mutate("ALLELE" = paste0(ALLELE," (",N_IND,")")) %>% #adding the number of possessors to allele labels
      mutate("ALLELE" = factor(ALLELE,levels = unique(.$ALLELE))) #changing allele variable back to factor
    
    p_list[[length(p_list)+1]] <- ggplot(dt)+
      geom_tile(aes(x = ALLELE,y = ID,fill = IS_IN_MOTHER,alpha = MARKER),col = "black")+
      theme_bw()+
      scale_fill_manual(values = c("TRUE" = "darkgreen","FALSE" = "red4"))+
      scale_alpha_manual(values = c("Iex2" = 1,"Iex3" = 0.6))+
      theme(axis.text.x = element_markdown(angle = 45,hjust = 1),panel.grid = element_blank())+
      labs(title = paste("MOTHER",i,"-",paste(marker_var,collapse = "_")), y = "OFFSPRING ID", 
           x = "ALLELE (POSSESSORS IN WHOLE DATASET)",fill = "IS IN MOTHER?")
  }
  p_final <- ggpubr::ggarrange(p_list[[1]],p_list[[2]],ncol = 1)
  
  ggsave(filename = paste0("./OUTCOMES/MY_PLOTS/after_filtering/mother_",i,"_",paste(marker_var,collapse = "_"),".png")
         ,device = "png",plot = p_final,width = 20,height = 14)
}
#
### SOLVING MOTHER'S HAPLOTYPES ####
#dataset
dt <- filtered_binary_dt %>% filter(PEDIGREE == "offspring")

#looping
out_list <- list()
for(i in sort(unique(dt$ID_MOTHER))){ #for each mother
  temp <- filter(dt,ID_MOTHER == i & IS_IN_MOTHER == T) %>% select(ID,ALLELE) %>% mutate("GENOTYPE" = 1) #subset for each mother
  #distance matrix
  distance_mat <- pivot_wider(temp,names_from = ALLELE,values_from = GENOTYPE,values_fill = 0) %>% arrange(ID)
  ids <- distance_mat$ID
  distance_mat <- distance_mat %>% select(-ID) %>% dist(method = "binary") %>% as.matrix()
  colnames(distance_mat) <- ids; rownames(distance_mat) <- ids
  #clustering
  cl <- cluster::pam(distance_mat,2,diss = T)
  for(j in c(1,2)){
    #extracting clusters
    cl_1 <- names(cl$clustering[cl$clustering == j])
    temp_1 <- select(temp,-GENOTYPE) %>% filter(ID %in% cl_1) %>% group_by(ID) %>%
      summarise("N_ALLELES" = n_distinct(ALLELE),"SET" = paste(sort(ALLELE),collapse = " "))
    sets <- group_by(temp_1,SET) %>% summarise("N_IND" = n(),"FREQ" = round(n()/nrow(.),2))
    temp_1 <- left_join(temp_1,sets,by = "SET")
    
    #intersection between individuals within a cluster + tolerance of a single allele if absent in only one individual 
    test_list <- strsplit(temp_1$SET,split = " ")
    my_hap <- get_union_haplotype(test_list)
    hap <- data.frame("ID_MOTHER" = rep(i,times = length(my_hap)),"ALLELE" = my_hap,"HAPLOTYPE" = j)
    
    out_list[[length(out_list)+1]] <- hap
  }
}
out_list <- as.data.frame(do.call(rbind,out_list))
write.csv(out_list,"generated_data/Mothers_haplotypes_intersect.csv",row.names = F,fileEncoding = "UTF-8")
#
### UNIQUE MATERNAL HAPLOTYPES ####
mother_haps <- read.csv("generated_data/Mothers_haplotypes_intersect.csv") %>% arrange(ID_MOTHER,HAPLOTYPE,ALLELE) %>%
  group_by(ID_MOTHER,HAPLOTYPE) %>% summarise("N_ALLELLES" = n_distinct(ALLELE),"ALLELES" = paste(ALLELE,collapse = " "))
length(unique(mother_haps$ALLELES))
#

### PLOTTING PER MOTHER - FILTERED WITH MATERNAL HAPLOTYPES ####
marker_var <- c("Iex2","Iex3")

mother_haps <- read.csv("generated_data/Mothers_haplotypes_intersect.csv") %>% group_by(ID_MOTHER,ALLELE) %>%
  mutate("OCCURRENCE" = n()) %>% ungroup() %>%
  mutate("HAPLOTYPE" = case_when(OCCURRENCE == 2 ~ "both",T ~ as.character(HAPLOTYPE))) %>%
  select(-OCCURRENCE) %>% distinct()

dt_org <- filtered_binary_dt %>% 
  filter(PEDIGREE == "offspring" & MARKER %in% marker_var) %>%
  select(ID,ALLELE,ID_MOTHER,IS_IN_MOTHER,MARKER) %>%
  group_by(ALLELE) %>% mutate("N_IND" = n_distinct(ID)) %>% #number of all individuals possessing a given allele
  left_join(mother_haps,by = c("ID_MOTHER","ALLELE")) %>% #joining mothers haplotypes
  mutate("HAPLOTYPE" = case_when(is.na(HAPLOTYPE) ~ "none",T ~ HAPLOTYPE)) #give label to paternal alleles

for(i in unique(dt_org$ID_MOTHER)){ #for all mothers
  p_list <- list() #list for plots
  dt_m <- dt_org %>% filter(ID_MOTHER == i)
  
  for(j in c("M","F")){ #grouped by mother or father
    dt <- dt_m %>% mutate("ID" = factor(as.character(ID),levels = clustered_order_m_or_f(.,parent = j))) #sorts individuals depending on maternal or paternal alleles similarity
    
    dt <- dt %>% arrange(desc(IS_IN_MOTHER),MARKER,desc(N_IND)) %>% #sorting alleles: maternal first then most frequent
      mutate("ALLELE" = paste0(ALLELE," (",N_IND,")")) %>% #adding the number of possessors to allele labels
      mutate("ALLELE" = case_when(HAPLOTYPE == "both" ~ paste0("<i style='color:blue4'>",ALLELE,"</i>"),
                                  HAPLOTYPE == "1" ~ paste0("<i style='color:green'>",ALLELE,"</i>"),
                                  HAPLOTYPE == "2" ~ paste0("<i style='color:red'>",ALLELE,"</i>"),
                                  HAPLOTYPE == "none" ~ paste0("<i style='color:black'>",ALLELE,"</i>"))) %>% #adding color to x-labels
      mutate("ALLELE" = factor(ALLELE,levels = unique(.$ALLELE))) #changing allele variable back to factor
    
    p_list[[length(p_list)+1]] <- ggplot(dt)+
      geom_tile(aes(x = ALLELE,y = ID,fill = IS_IN_MOTHER,alpha = MARKER),col = "black")+
      theme_bw()+
      scale_fill_manual(values = c("TRUE" = "darkgreen","FALSE" = "red4"))+
      scale_alpha_manual(values = c("Iex2" = 1,"Iex3" = 0.6))+
      theme(axis.text.x = element_markdown(angle = 45,hjust = 1),panel.grid = element_blank())+
      labs(title = paste("MOTHER",i,"-",paste(marker_var,collapse = "_")), y = "OFFSPRING ID", 
           x = "ALLELE (POSSESSORS IN WHOLE DATASET)",fill = "IS IN MOTHER?")
  }
  p_final <- ggpubr::ggarrange(p_list[[1]],p_list[[2]],ncol = 1)
  
  ggsave(filename = paste0("./OUTCOMES/MY_PLOTS/filtered_and_maternalHaplotypes/mother_",i,"_",paste(marker_var,collapse = "_"),".png")
         ,device = "png",plot = p_final,width = 20,height = 14)
}
#