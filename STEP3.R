### HAPLOTYPE RECONSTRUCTION ###
### Tomasz Gaczorek ###
### tomasz.gaczorek@doctoral.uj.edu.pl ###
### 11.01.2023 ###

### INTRO ####
#STEP 3 - reconstruction of paternal haplotypes
#REQUIRE: Mothers_haplotypes_intersect.csv (generated in STEP 1); Lissotriton_ids_locality.csv

#setwd("~/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG/annealing/RDSs/")
setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG")

library(dplyr)
#
### FUNCTIONS ####
find_similar <- function(vec1,list_v){
  x <- which.min(sapply(list_v, function(x){sum(!(vec1 %in% x))+sum(!(x %in% vec1))}))
  names(x)[1]
} #find the list element which is the most similar to provided vector

min_max_hap <- function(list_v){
  lens <- sapply(list_v,length)
  min_ind <- which.min(lens)
  max_ind <- which.max(lens)
  list(list_v[[min_ind]],list_v[[max_ind]])
} #find the shortest and the longest element of a list

min_max_comb_hap <- function(list_v){
  lens <- sapply(list_v,length)
  min_ind <- which.min(lens); min_value <- lens[min_ind]
  max_ind <- which.max(lens); max_value <- lens[max_ind]
  
  if(min_value != max_value){
    diff <- max_value-min_value
    n_comb <- sum(sapply(c(1:diff),choose,n = diff))+1
    return(list("MIN" = list_v[[min_ind]],"MAX" = list_v[[max_ind]],"N_COMB" = n_comb))
  } else {
    return(list("MIN" = list_v[[min_ind]],"MAX" = NA,"N_COMB" = 1))
  }
} #length of the shortest and the longest element of list + calculation of total number of possible combinations

get_mother_hap_rep_ll <- function(x){
  splited <- strsplit(x,"_|\\.")[[1]]
  ll <- readRDS(x)$fit
  return(c(splited[c(2,6,4)],ll))
} #gathering data from RDS files

does_contain_more <- function(list_v,shorter){ 
  logicals <- sapply(list_v,function(x){sum(shorter %in% x) == length(shorter)})
  logicals
} #checking which list elements contain shorter vector
#
### DATASETS ####
df_mother_hap <- read.csv("generated_data/Mothers_haplotypes_intersect.csv")

mothers <- unique(df_mother_hap$ID_MOTHER)
#

### COLLECTING FROM RDS ####
rds <- list.files("OUTCOMES/RDSs/",pattern = "\\.")[list.files(".",pattern = "\\.") != "log.txt"]
rds <- paste0("OUTCOMES/RDSs/",rds)
dir.create("OUTCOMES/Results")

#gathering LL values
lls <- as.data.frame(t(sapply(rds,get_mother_hap_rep_ll)))
rownames(lls) <- NULL; colnames(lls) <- c("MOTHER","HAP_NUM","REP","LL")
lls <- arrange(lls,MOTHER,HAP_NUM,REP) %>% mutate("LL" = as.numeric(LL))
write.csv(lls,"OUTCOMES/Results/lls.csv",row.names = F,fileEncoding = "UTF-8")

#lls <- read.csv("lls.csv") #backdoor

#collecting relevant results into a single list
eternal_list <- list()
for(j in mothers){
  lls_temp <- filter(lls,MOTHER == j) %>%
    filter(LL >= max(LL)-2) %>%
    filter(HAP_NUM == min(HAP_NUM)) #choice of the best number of haplotypes
  
  useful_rds <- paste0("OUTCOMES/RDSs/MOTHER_",j,"_REP_",lls_temp$REP,"_HAP_",unique(lls_temp$HAP_NUM),".RDS") #proper RDS
  
  #gathering alternatives from the equally good replicates
  alternatives <- list()
  for(i in useful_rds){
    alts_temp <- readRDS(i)$alternative_conf
    alternatives[(length(alternatives)+1):(length(alternatives)+length(alts_temp))] <- alts_temp
  }
  
  #template
  template <- alternatives[[1]]
  names(template) <- paste0("HAP_",c(1:length(template)))

  #giving names & sorting
  for(i in 1:length(alternatives)){
    chosen_names <- sapply(alternatives[[i]],find_similar,template)
    names(alternatives[[i]]) <- chosen_names
    alternatives[[i]] <- alternatives[[i]][sort(chosen_names)]
  }
  
  #taking out corresponding haplotypes
  haps <- list()
  for(i in 1:length(template)){
    haps[[i]] <- lapply(alternatives,"[[",i)
  }

  #filtering duplicated alternatives
  haps_filtered <- lapply(haps, unique)

  #return shortest, longest and number of combinations
  p <- lapply(haps_filtered,min_max_comb_hap)
  p <- lapply(p,function(x){x$MOTHER <- j;x})
  
  #gathering for each mother
  eternal_list[[as.character(j)]] <- p
}
saveRDS(eternal_list,"OUTCOMES/Results/Results.RDS")

#creating summary
out_list <- list()
for(i in 1:length(eternal_list)){
  for(j in 1:length(eternal_list[[i]])){
    out_list[[length(out_list)+1]] <- c(names(eternal_list)[i],length(eternal_list[[i]]),eternal_list[[i]][[j]]$N_COMB)
  }
}
out_list <- as.data.frame(do.call(rbind,out_list))
colnames(out_list) <- c("MOTHER","N_HAPS","N_COMB")
write.csv(out_list,"OUTCOMES/Results/Summary.csv",row.names = F,fileEncoding = "UTF-8")

### ANALYSIS OF COLLECTED DATA  - HAS TO BE MODIFIED DEPENDING ON WHAT WE ARE GOING TO DO####
results <- readRDS("OUTCOMES/Results/Results.RDS")

#group paternal haplotypes from multiple mothers
#establish mothers
m_results <- names(results)
haps_lens <- sapply(results,length) 
mothers <- unlist(sapply(c(1:length(m_results)), function(x){rep(m_results[x],times = haps_lens[x])}))

#checking overlapping
shortest <- unlist(lapply(results,function(y){lapply(y,"[[","MIN")}),recursive = F,use.names = F)
longest <- unlist(lapply(results,function(y){lapply(y,"[[","MAX")}),recursive = F,use.names = F)
longest[is.na(longest)] <- shortest[is.na(longest)]

x <- shortest

my_list<- list()
while(length(x)>0){ #grouping paternal haplotypes
  repeated <- does_contain_more(longest,x[[1]])
  my_list[[length(my_list)+1]] <- list("MIN"=x[repeated],"MAX" = longest[repeated],"MOTHERS" = mothers[repeated])
  #indexes <- which(repeated == F)
  x <- x[!repeated]
  longest <- longest[!repeated]
  mothers <- mothers[!repeated]
}

saveRDS(my_list,"OUTCOMES/Results/Groups_paternalHaps.RDS")

# gather info about presumable populations
groups_hap <- lapply(c(1:length(my_list)),function(x){data.frame("ID" = rep(x,times = length(my_list[[x]]$MOTHERS)),
                                                                 "MOTHER" = my_list[[x]]$MOTHERS)
  }) 
groups_hap <- as.data.frame(do.call(rbind,groups_hap))

localities <- read.csv("generated_data/Lissotriton_ids_locality.csv") %>% select(INDIVIDUAL_ID,LOCALITY,SPECIES) %>%
  mutate("INDIVIDUAL_ID" = as.character(INDIVIDUAL_ID))

groups_hap <- left_join(groups_hap,localities,by = c("MOTHER"="INDIVIDUAL_ID"))
#write.csv(groups_hap,"HapGroups_Localities.csv",row.names = F,fileEncoding = "UTF-8")

#localities mapping
map <- read.csv("generated_data/Lissotriton_ids_locality.csv") %>%
  filter(LOCALITY %in% unique(groups_hap$LOCALITY)) %>% select(LOCALITY,LATITUDE,LONGITUDE) %>% distinct()
write.csv(map,"3ponds_map.csv",row.names = F,fileEncoding = "UTF-8")
  
#
