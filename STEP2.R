### HAPLOTYPE RECONSTRUCTION ###
### Tomasz Gaczorek ###
### tomasz.gaczorek@doctoral.uj.edu.pl ###
### 11.01.2023 ###

### INTRO ####
#STEP 2 - reconstruction of paternal haplotypes
#REQUIRE: binary_dt_filtered.csv; Mothers_haplotypes_intersect.csv (generated in STEP 1)

#("~/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG/annealing/")
setwd("C:/Users/230982/Dropbox/MHC_projekt/Shared_MHC_project/Lissotriton_MHC_haplotypes_R/TG")
time_start <- Sys.time() #checking script duration
library(dplyr)
library(doParallel)
library(foreach)
#
### FUNCTIONS ####
which.haplotype <- function(offspring_tab,hap_table,child){ #checks which of maternal haplotypes you possess
  mother_dt <- filter(hap_table,ID_MOTHER == unique(offspring_tab$ID_MOTHER[offspring_tab$ID == child]))
  matches <- c(sum(offspring_tab$ALLELE[offspring_tab$ID == child] %in% mother_dt$ALLELE[mother_dt$HAPLOTYPE == 1]),
               sum(offspring_tab$ALLELE[offspring_tab$ID == child] %in% mother_dt$ALLELE[mother_dt$HAPLOTYPE == 2]))
  c(1,2)[which.max(matches)]
} #finding which maternal haplotype an offspring has

likelihood <- function(child_alleles,mother_alleles,compared_father,error){
    in_mother <- child_alleles %in% mother_alleles
    in_father <- child_alleles %in% compared_father
    only_father <- compared_father[!(compared_father %in% child_alleles)]
    only_mother <- mother_alleles[!(mother_alleles %in% child_alleles)]
    out <- rep(1-error,times = length(child_alleles)) # 1-a unless
    out[!in_mother & !in_father] <- error # absent in both mother and father
    #out[in_mother & in_father] <- 1-error #present in both
    out <- c(out,rep(error,times = length(only_father)),rep(error,times = length(only_mother))) #add error for every paternal and expected maternal allele absent in a child
    return(prod(out)) 
} #calculate likelihood across all offspring's alleles

get_list_childHap_fromMother <- function(vec_which_hap,df_mother_hap,df_offspring_tab){
  out_list <- list()
  for(i in names(vec_which_hap)){
    mother_no <- unique(df_offspring_tab$ID_MOTHER[df_offspring_tab$ID == i])
    out_list[[i]] <- df_mother_hap$ALLELE[df_mother_hap$HAPLOTYPE == vec_which_hap[i] & df_mother_hap$ID_MOTHER == mother_no]
  }
  out_list
} #list of maternal haplotypes for the offspring

overall_likelihood <- function(conf,list_child_haps,list_child_fromMother,error){
  sum_likeli <- numeric()
  for(i in names(list_child_haps)){
    sum_likeli[length(sum_likeli)+1] <- max(sapply(conf,likelihood,child_alleles = list_child_haps[[i]],
                                                   mother_alleles = list_child_fromMother[[i]],error = error)) #likelihood of configuration for one individual
  }
  sum(log(sum_likeli))
} #likelihood across many children

random_list_elements <- function(my_list,size){
  all_elements_lengths <- sapply(my_list,length)
  cum_length <- cumsum(all_elements_lengths)
  chosen <- sample(c(1:sum(all_elements_lengths)),size = size,replace = F)
  elements <- sapply(chosen, function(x){min(which(cum_length >= x))})
  inside <- sapply(chosen, function(x){temp <-cum_length-x;
  if(sum(temp == 0) != 0){ #if the last element of a vector
    out <- all_elements_lengths[which(temp == 0)]
  } else if(sum(temp < 0) == 0){ #if coming from the first vector
    out <- temp[1]
  } else {
    out <- abs(max(temp[temp<0])) #otherwise
  }
  return(out)})
  return(list("vectors" = elements,"coordinates" = inside))
} #randomly selects element from the list

delete_alleles <- function(conf,size,minimum_length = 3){
  temp <- conf
  if(sum(sapply(conf,length)) != 0){
    indexes <- random_list_elements(my_list = conf,size = size)
    for(i in 1:length(indexes$vectors)){
      if(length(temp[[indexes$vectors[i]]]) > minimum_length){
        temp[[indexes$vectors[i]]] <- temp[[indexes$vectors[i]]][-(indexes$coordinates[i])] 
      }
    }
  }
  temp
} #deletes alleles from the list

add_alleles <- function(conf,size,set_alleles){
  which_haplotypes <- sample(c(1:length(conf)),size = size,replace = T)
  temp <- conf
  for(i in which_haplotypes){
    if(length(set_alleles[!(set_alleles %in% temp[[i]])]) > 0){
      temp[[i]] <- c(temp[[i]],sample(set_alleles[!(set_alleles %in% temp[[i]])],size = 1))
    }
  }
  temp
} #add alleles to the list

fight_with_NULLs <- function(conf){
  temp <- conf
  if(length(conf)>1){
    nullek <- which(sapply(conf,is.null))[1]
    x <- c(1:length(conf))
    from <- sample(x[x != nullek],1)
    if(is.null(temp[[from]]) == F){
      temp[[nullek]] <- temp[[from]]
    }
  }
  temp
} #replaces empty haplotypes with one randomly selected

running <- function(init_conf, iter=20000,offspring_tab,mother_tab,error,n_hap){
  
  count <- 1
  #initialization of explored solution and best solution
  temp_conf <- init_conf
  best_conf <- init_conf
  
  #specifying haplotypes
  alleles_set <- offspring_tab %>% pull(ALLELE) %>% unique() 
  child_haplotypes <- split(f = offspring_tab$ID,x = offspring_tab$ALLELE) #children haplotypes [named list]
  which_mother_hap <- sapply(names(child_haplotypes),which.haplotype,offspring_tab = offspring_tab,hap_table = mother_tab) #which mother's haplotype the child possess [named vector]
  child_haplotypes_fromMother <- get_list_childHap_fromMother(vec_which_hap = which_mother_hap,
                                                              df_mother_hap = mother_tab,
                                                              df_offspring_tab = offspring_tab) #children's haplotypes of the maternal origin [named list]
  #setting fit values
  fit <- overall_likelihood(conf = temp_conf,list_child_haps = child_haplotypes,
                            list_child_fromMother = child_haplotypes_fromMother,error = error)
  bestfit <- fit
  
  #sinking
  alternative_confs <- list(best_conf)
  l_hoods <- bestfit
  
  ## the simulated annealing loop
  count <- 1
  while(count < iter){
    
    test_conf <- temp_conf
    if(count > (iter/2)){
      if(sum(sapply(test_conf,is.null))>0){
        test_conf <- fight_with_NULLs(test_conf) #avoiding an algorithm from being stuck in case of paternal haplotypes similarity
      }
    }
    
    #changes
    if(count <= (.8*iter)){
      choice <- sample(c("D","A"),size = 1,prob = c(.5,.5))
    } else {
      choice <- "A"
    }
    
    if(choice == "D"){
      test_conf <- delete_alleles(test_conf,1)
    } else if(choice == "A"){
      test_conf <- add_alleles(test_conf,1,set_alleles = alleles_set)
    } 
    
    #obtaining the testing solution x'
    if(identical(test_conf,temp_conf) == F){
      testfit <- overall_likelihood(conf = test_conf,list_child_haps = child_haplotypes,
                                    list_child_fromMother = child_haplotypes_fromMother,error = error)
    } else {
      testfit <- fit
    }
    
    #checking if we replace x by x'
    if((fit-testfit) < 0){ #always accept if better 
      temp_conf <- test_conf
      fit <- testfit
    } else if(((fit-testfit) == 0) & (count > (.8*iter))){ #accept if equally good after 80% of cycles
      temp_conf <- test_conf
      fit <- testfit
    }
    
    #updating the best solution
    if(testfit > bestfit){ #if better
      
      #sinking
      if(sum(sapply(alternative_confs,identical,y = test_conf))==0){ #if an alternative hasn't been collected yet
        indexes_in_2LL <- which((testfit-l_hoods)<=2)
        l_hoods <- c(l_hoods[indexes_in_2LL],testfit)
        alternative_confs <- alternative_confs[indexes_in_2LL]
        alternative_confs[[length(alternative_confs)+1]] <- lapply(test_conf,sort) #collect an alternative and any previously collected alternatives within 2LL range 
      }
      best_conf <- test_conf
      bestfit <- testfit
      count <- 1 #start from beginning
    } else if(testfit == bestfit){ #if equally good 
      count <- count + 1 #extend a chain
      if(sum(sapply(alternative_confs,identical,y = test_conf))==0){
        l_hoods <- c(l_hoods,testfit)
        alternative_confs[[length(alternative_confs)+1]] <- lapply(test_conf,sort) #collect an alternative if not present before
      }
    } else{ #if worse
      count <- count + 1 #extend a chain
    }
  }
  #returning the solution
  return(list("conf"=lapply(best_conf,sort), "fit"=bestfit,"L_hoods" = l_hoods,"alternative_conf" = alternative_confs))
} #algorithm for finding optimal paternal configuration
#
### DATASETS ####
df_offspring <- read.csv("generated_data/binary_dt_filtered.csv",stringsAsFactors = F) %>% #data with alleles
  filter(PEDIGREE == "offspring" & MARKER %in% c("Iex2","Iex3")) %>%
  select(ID,ALLELE,ID_MOTHER)

df_mother_hap <- read.csv("generated_data/Mothers_haplotypes_intersect.csv",stringsAsFactors = F) #data with maternal haplotypes

mothers <- unique(df_mother_hap$ID_MOTHER) #mothers' IDs
#
### BODY ####
#constant
replicates <- 1:3 

#setting up paralellization
no_cores <- 10
cl <- makeCluster(no_cores)
registerDoParallel(cl)

dir.create("OUTCOMES/RDSs")
writeLines("","OUTCOMES/RDSs/log.txt")

invisible(foreach(m = mothers,.packages = "dplyr") %:%
            foreach(k = replicates,.packages = "dplyr") %dopar% { #all mothers and replicates in parallel
              
              df_offspring_v <- filter(df_offspring,ID_MOTHER == m)
              df_mother_hap_v <- df_mother_hap %>% filter(ID_MOTHER == m)
              
              num_hap_estimated <- filter(df_offspring_v,!(ALLELE %in% df_mother_hap_v$ALLELE)) %>% select(-ID_MOTHER) %>%
                group_by(ID) %>% summarise("SET" = paste(sort(ALLELE),collapse = ",")) %>% distinct(SET) %>% nrow() #number of unique paternal haplotypes
              
              haplotypes_number <- c((num_hap_estimated-3):(num_hap_estimated+1)) #set of numbers of paternal haplotypes 
              haplotypes_number <- haplotypes_number[haplotypes_number > 0] #must be larger than 0
              
              for(i in haplotypes_number){
                #running
                initial_configuration <- replicate(i,c(NULL)) #empty haplotypes
                temp <- running(init_conf = initial_configuration,offspring_tab = df_offspring_v,
                                mother_tab = df_mother_hap_v,error = 0.03,n_hap = i,iter = i*20000) #running searching algorithm
                
                if(k %% 1 == 0){ #in every step
                  sink("OUTCOMES/RDSs/log.txt",append = T)
                  print(paste("MOTHER",m,"REP",k,"HAP_NUM",i)) #update log
                  sink()
                  rds_name <- paste("OUTCOMES/RDSs/MOTHER_",m,"_REP_",k,"_HAP_",i,".RDS",sep = "") #save outcome to RDS file
                  saveRDS(temp,file = rds_name)
                }
              }
            })

stopCluster(cl)
print(Sys.time()-time_start)
#

