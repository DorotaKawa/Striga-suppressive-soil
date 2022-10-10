setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/POP#2/GJAM/2021Dec-redo/")

library(readxl)
library(openxlsx)
library(tidyverse)
library(naniar)
library(pheatmap)

path<- "20220625-Residual Correlation_M_complete.xlsx"
tab_names <- excel_sheets(path = path)
list_all <- lapply(tab_names, function(x) read_excel(path = path, sheet = x, na= 'NA'))
lapply(list_all, head, n=2)
names(list_all) <- tab_names
names(list_all) <- gsub("output", "", names(list_all))

#### Filter meaningful (=logical for each phenotype) corelations (different traits for 2 and 3 weeks, so calculated separately)

FunCor2 <- list_all[grep("Fun.*2week", names(list_all))]
lapply(FunCor2, head, n=2)

FunCor3 <- list_all[grep("Fun.*3week", names(list_all))]
lapply(FunCor3, head, n=2)

#### Fungi 2week - filter positive/negative correlations per trait
selcor2 <- function(x) {
  x %>% 
    replace_with_na_at(.vars = c("DMBQ", "Strigaattachmentspergroot", "Syringicacid", "Vanillicacid"),
                       condition = ~.x > 0) %>%
    replace_with_na_at(.vars = c("Aerenchyma"),
                       condition = ~.x < 0)
}

Fun2_list <- FunCor2 %>% lapply(selcor2)

##### Fungi 3week - filter positive/negative correlations per trait
selcor3 <- function(x) {
  x %>% 
    replace_with_na_at(.vars = c("Strigaattachmentspergroot"),
                       condition = ~.x > 0) %>%
    replace_with_na_at(.vars = c("Aerenchyma", "Suberin"),
                       condition = ~.x < 0)
}

Fun3_list <- FunCor3 %>% lapply(selcor3)
Fun_list <- append(Fun2_list, Fun3_list)

#################################### First make ranking for individual trait across all microbial compartments ##################################################

#### Create one df for each trait combinig values for each microbial compartment

contains_tof <- function(x, nam) sum(grepl(nam, colnames(x))) == 1
sel_names <- c("term", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

### For 2 week traits 
tof2 <- c("Aerenchyma", "DMBQ", "Strigaattachmentspergroot", "Syringicacid", "Vanillicacid")  
tofL2 <- vector(length = length(tof2), mode='list')
names(tofL2) <- tof2
for(i in tof2){
  Fun <- Fun_list[sapply(Fun_list, contains_tof, nam = i)]
  t.sel_names <- c(sel_names, i)
  trait <- lapply(Fun, "[", t.sel_names)
  
  trait2 <- lapply(1:length(trait), function(x) { 
    tL <- trait [[x]]
    colnames(tL) <- c("term","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", names(trait)[x])
    return(tL)})
  
  tofL2[[i]] <- Reduce(
    function(x, y, ...) full_join(x, y, by = c('term',"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
    trait2)
}
lapply(tofL2, head, n=2)

# write.xlsx(tofL2, "Fungi - Correlations 2 weeks per trait.xlsx")

### For 3 week traits 
tof3 <- c("Aerenchyma", "Strigaattachmentspergroot", "Suberin")  
tofL3 <- vector(length = length(tof3), mode='list')
names(tofL3) <- tof3
for(i in tof3){
  Fun <- Fun_list[sapply(Fun_list, contains_tof, nam = i)]
  t.sel_names <- c(sel_names, i)
  trait <- lapply(Fun, "[", t.sel_names)
  
  trait2 <- lapply(1:length(trait), function(x) { 
    tL <- trait [[x]]
    colnames(tL) <- c("term","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", names(trait)[x])
    return(tL)})
  
  tofL3[[i]] <- Reduce(
    function(x, y, ...) full_join(x, y, by = c('term',"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
    trait2)
}
lapply(tofL3, head, n=2)

# write.xlsx(tofL3, "Fungi - Correlations 3 weeks per trait.xlsx")

### Collect values from each compartment per trait 

Comb_comps <- function(x) {
  x %>%
    gather(key = "comp", value = "corr", -(term:Species))
}

tofL2comb <- tofL2 %>% lapply(Comb_comps)
lapply(tofL2comb, head, n=2)

tofL3comb <- tofL3 %>% lapply(Comb_comps)
lapply(tofL3comb, head, n=2)

tofL2comb_rename <- list()
for(i in 1:length(tofL2comb)){
  tname    <- names(tofL2comb)[i]
  tdf      <- tofL2comb[[i]]
  newnames <- rep('', ncol(tdf))
  newnames[length(newnames)] <- paste0('_', tname)
  colnames(tdf) <- paste0(colnames(tdf), newnames)
  tofL2comb_rename[[i]] <- tdf
}
names(tofL2comb_rename) <- names(tofL2comb)
lapply(tofL2comb_rename, head, n=3)


tofL3comb_rename <- list()
for(i in 1:length(tofL3comb)){
  tname    <- names(tofL3comb)[i]
  tdf      <- tofL3comb[[i]]
  newnames <- rep('', ncol(tdf))
  newnames[length(newnames)] <- paste0('_', tname)
  colnames(tdf) <- paste0(colnames(tdf), newnames)
  tofL3comb_rename[[i]] <- tdf
}
names(tofL3comb_rename) <- names(tofL3comb)
lapply(tofL3comb_rename, head, n=3)

#### Merge all df per timepoint into one
#Corr2 and Corr3 have correlation values for each trait combined across compartments

Corr2 <- Reduce(
  function(x, y, ...) full_join(x, y, by = c('term',"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "comp")),
  tofL2comb_rename)

Corr3 <- Reduce(
  function(x, y, ...) full_join(x, y, by = c('term',"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "comp")),
  tofL3comb_rename)

Corr23 <- full_join(select(Corr2, term:comp,corr_Aerenchyma, corr_Strigaattachmentspergroot),
                    select(Corr3, term:comp,corr_Aerenchyma, corr_Strigaattachmentspergroot))
                      
###Create individual ranking for each trait across compartments

Rank_allcomp2 <- Corr2 %>% 
  mutate('DMBQ_rank' = rank(corr_DMBQ, na.last = "keep")) %>%
  mutate('Syringic_rank' = rank(corr_Syringicacid, na.last = "keep")) %>%
  mutate('Vanillic_rank' = rank(corr_Vanillicacid, na.last = "keep")) 

Rank_allcomp3 <- Corr3 %>% 
mutate('Sub_rank' = rank(-corr_Suberin, na.last = "keep")) 

Rank_allcomp23 <- Corr23 %>%
  mutate('Aer_rank' = rank(-corr_Aerenchyma, na.last = "keep")) %>%
  mutate("Att_rank" = rank(corr_Strigaattachmentspergroot, na.last = "keep"))

### Since less fungal taxa found, instead of top 100, focus on top 50
top_2w <- Rank_allcomp2 %>%
  replace_with_na_at(vars(DMBQ_rank:Vanillic_rank), condition = ~.x > 50) %>%
  filter_at(vars(-(term:corr_Vanillicacid)), any_vars(! is.na(.))) 

top_3w <- Rank_allcomp3 %>%
  replace_with_na_at(vars(Sub_rank), condition = ~.x > 50) %>%
  filter_at(vars(-(term:corr_Suberin)), any_vars(! is.na(.))) 

top_23w <- Rank_allcomp23 %>%
  replace_with_na_at(vars(Aer_rank:Att_rank), condition = ~.x > 50) %>%
  filter_at(vars(-(term:corr_Strigaattachmentspergroot)), any_vars(! is.na(.))) 

### 2wpi
tof2 <- c("DMBQ_rank", "Syringic_rank", "Vanillic_rank")
hmsel2List <- vector(length = length(tof2), mode = "list")
names(hmsel2List) <- tof2

for (i in tof2) {
  hmsel2List[[i]] <- top_2w %>%
    select(term:comp, i) %>%
    spread(comp, i) %>%
    select(term:Species, contains("2week")) %>%
    arrange(Phylum) %>%
    mutate(sum = rowSums(dplyr::select(., FunBK2week, FunEndo2week, FunRh2week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum) %>%
    select(term:Species,FunBK2week, FunRh2week, FunEndo2week)
  
}

pdf("Fungi top 50 heatmap 2wpi.pdf")
for (i in tof2){
  hp <- hmsel2List[[i]]
  hp_m <- data.matrix(hp[,9:11])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 5, color = mycolors, show_rownames = FALSE,
             annotation_row=annotation, main =i )
  )
}
dev.off()

### 3wpi
tof3 <- c("Sub_rank")
hmsel3List <- vector(length = length(tof3), mode = "list")
names(hmsel3List) <- tof3

for (i in tof3) {
  hmsel3List[[i]] <- top_3w %>%
    select(term:comp, i) %>%
    spread(comp, i) %>%
    select(term:Species, contains("3week")) %>%
    arrange(Phylum) %>%
    mutate(sum = rowSums(dplyr::select(., FunBK3week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum)
  }

pdf("Fungi top 50 heatmap 3wpi.pdf")
for (i in tof3){
  hp <- hmsel3List[[i]]
  hp_m <- data.matrix(hp[,9])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 5, color = mycolors, show_rownames = FALSE,
             annotation_row=annotation, main =i )
  )
}
dev.off()

### 2 and 3wpi
tof23 <- c("Att_rank", "Aer_rank")
hmsel23List <- vector(length = length(tof23), mode = "list")
names(hmsel23List) <- tof23

for (i in tof23) {
  hmsel23List[[i]] <- top_23w %>%
    select(term:comp, i) %>%
    spread(comp, i) %>%
    select(term:Species, contains("week")) %>%
    arrange(Phylum) %>%
    mutate(sum = rowSums(dplyr::select(.,FunBK2week, FunBK3week, FunEndo3week, FunEndoSand3week, FunRh3week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum) %>%
    select(term:Species,FunBK2week, FunBK3week, FunRh3week, FunEndo3week, FunEndoSand3week)
}


pdf("Fungi top 50 heatmap 2 and 3wpi.pdf")
for (i in tof23){
  hp <- hmsel23List[[i]]
  hp_m <- data.matrix(hp[,9:13])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 5, color = mycolors, show_rownames = FALSE,
             annotation_row=annotation, main =i )
  )
}
dev.off()


