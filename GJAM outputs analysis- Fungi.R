setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/POP#2/GJAM/2021Dec-redo/")

library(readxl)
library(openxlsx)
library(tidyverse)
library(naniar)
library(pheatmap)
library(wesanderson)
mycolors <- wes_palette("Zissou1", 10, type = "continuous")

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

#################################### First make rankings for individual trait across all microbial compartments ##################################################

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

### Select top 50 taxa to plot
top_2w <- Rank_allcomp2 %>%
  replace_with_na_at(vars(DMBQ_rank:Vanillic_rank), condition = ~.x > 50) 

top_3w <- Rank_allcomp3 %>%
  replace_with_na_at(vars(Sub_rank), condition = ~.x > 50)  

top_23w <- Rank_allcomp23 %>%
  replace_with_na_at(vars(Aer_rank:Att_rank), condition = ~.x > 50) 

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
    mutate(sum = rowSums(dplyr::select(., FunBK2week, FunRh2week, FunEndo2week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum) %>%
    select(term:Species, FunBK2week, FunRh2week, FunEndo2week)
  
}

pdf("Fun top 50 heatmap 2wpi.pdf")
for (i in tof2){
  hp <- hmsel2List[[i]]
  hp_m <- data.matrix(hp[,9:11])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 4, fontsize = 8, color = mycolors, show_rownames = FALSE,
             annotation_row=annotation, main =i )
  )
}
dev.off()

### 3wpi
tof3 <- c("Sub_rank")
hmsel3List <- vector(length = length(tof3), mode = "list")
names(hmsel3List) <- tof3

for (i in tof3) {
  hmsel3List[[i]]  <- top_3w %>%
    select(term:comp, i) %>%
    spread(comp, i) %>%
    select(term:Species, contains("3week")) %>%
    arrange(Phylum) %>%
    mutate(sum = rowSums(dplyr::select(., FunBK3week, FunRh3week, FunEndo3week, FunEndoSand3week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum)%>%
    select(term:Species, FunBK3week, FunRh3week, FunEndo3week, FunEndoSand3week)
}

pdf("Fun top 50 heatmap 3wpi.pdf")
for (i in tof3){
  hp <- hmsel3List[[i]]
  hp_m <- data.matrix(hp[,9:12])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 4, fontsize = 8, color = mycolors, show_rownames = FALSE,
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
    mutate(sum = rowSums(dplyr::select(.,FunBK2week, FunRh2week, FunEndo2week, FunBK3week, FunRh3week, FunEndo3week, FunEndoSand3week), na.rm = TRUE)) %>%
    filter(sum != 0) %>%
    select(-sum) %>%
    select(term:Species,FunBK2week, FunRh2week, FunEndo2week,FunBK3week, FunRh3week, FunEndo3week, FunEndoSand3week)
}


pdf("Fun top 50 heatmap 2 and 3wpi.pdf")
for (i in tof23){
  hp <- hmsel23List[[i]]
  hp_m <- data.matrix(hp[,9:15])
  annotation <- data.frame(row.names =  hp$term, "Phylum" = hp$Phylum, "Class" = hp$Class)
  rownames(hp_m) <-hp$term
  
  print(
    pheatmap(hp_m, cluster_cols = FALSE, cluster_rows = FALSE, clustering_distance_rows = "euclidean", 
             border_color = "gray80", cellwidth = 10, cellheight = 4, fontsize = 8, color = mycolors, show_rownames = FALSE,
             annotation_row=annotation, main =i )
  )
}
dev.off()


#################################### Combined rannkings for individual traits across all microbial compartments ##################################################
##### Calculate individual and combined ranks for each compartment separately

ranking2 <- function(x) {
  x %>% 
    mutate('Aer_rank' = rank(-Aerenchyma, na.last = "keep")) %>%
    mutate('DMBQ_rank' = rank(DMBQ, na.last = "keep")) %>%
    mutate('Syringic_rank' = rank(Syringicacid, na.last = "keep")) %>%
    mutate('Vanillic_rank' = rank(Vanillicacid, na.last = "keep")) %>%
    mutate("Att_rank" = rank(Strigaattachmentspergroot, na.last = "keep")) %>%
    mutate("Aer_combined_rank" = rank((Aer_rank + Att_rank), na.last = "keep")) %>%
    mutate("DMBQ_combined_rank" = rank((DMBQ_rank + Att_rank), na.last = "keep")) %>%
    mutate("Syringic_combined_rank" = rank((Syringic_rank + Att_rank), na.last = "keep")) %>%
    mutate("Vanillic_combined_rank" = rank((Vanillic_rank + Att_rank), na.last = "keep")) %>%
    select(term:Species, contains('rank'))
}

ranking3 <- function(x) {
  x %>% 
    mutate('Aer_rank' = rank(-Aerenchyma, na.last = "keep")) %>%
    mutate('Sub_rank' = rank(-Suberin, na.last = "keep")) %>%
    mutate("Att_rank" = rank(Strigaattachmentspergroot, na.last = "keep")) %>%
    mutate("Aer_combined_rank" = rank((Aer_rank + Att_rank), na.last = "keep")) %>%
    mutate("Sub_combined_rank" = rank((Sub_rank + Att_rank), na.last = "keep"))%>%
    select(term:Species, contains('rank'))
}

Fun2_ranks <- Fun2_list %>% lapply (ranking2)
Fun3_ranks <- Fun3_list %>% lapply (ranking3)

Fun_ranks <- append(Fun2_ranks, Fun3_ranks)
# write.xlsx(Fun_ranks, "Ranks per trait*comp Fungi.xlsx")

##### Plot correlations in scatterplots and mark the combiranks values

#first get all residual correlations, not only the meaningful (=filtered) ones and combine them with the file with ranks

combo2 <- map2(FunCor2, Fun2_ranks, full_join)
lapply(combo2, head, n=2)

combo3 <- map2(FunCor3, Fun3_ranks, full_join)
lapply(combo3, head, n=2)


#scatter plots showing proportion of observations with combined rank 

### 3 week plots - suberin
pdf("3 Suberin combirank - Fungi.pdf")
for(i in names(combo3)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo3[[i]]), aes(Strigaattachmentspergroot,Suberin))+
      geom_point(aes(color = Sub_combined_rank)) +
      ggtitle(i)+
      xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  )
}
dev.off()

pdf("3 Aerenchyma combirank - Fungi.pdf")
for(i in names(combo3)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo3[[i]]), aes(Strigaattachmentspergroot,Aerenchyma))+
      geom_point(aes(color = Aer_combined_rank)) +
      ggtitle(i)+
      xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  )
}
dev.off()

### 2 week plots - aerenchyma and HIFs
pdf("2 Aerenchyma combirank - Fungi.pdf")
for(i in names(combo2)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo2[[i]]), aes(Strigaattachmentspergroot,Aerenchyma))+
      geom_point(aes(color = Aer_combined_rank)) +
      ggtitle(i)+
      xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  ) 
}
dev.off()

pdf("2 DMBQ combirank - Fungi.pdf")
for(i in names(combo2)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo2[[i]]), aes(Strigaattachmentspergroot,DMBQ))+
      geom_point(aes(color = DMBQ_combined_rank)) +
      ggtitle(i)+
      xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  )
}
dev.off()

pdf("2 VA combirank - Fungi.pdf")
for(i in names(combo2)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo2[[i]]), aes(Strigaattachmentspergroot,Vanillicacid))+
      geom_point(aes(color = Vanillic_combined_rank)) +
      ggtitle(i)+xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  )
}
dev.off()

pdf("2 SA combirank - Fungi.pdf")
for(i in names(combo2)){
  print(
    plotdf1 <-  ggplot(as.data.frame(combo2[[i]]), aes(Strigaattachmentspergroot,Syringicacid))+
      geom_point(aes(color = Syringic_combined_rank)) +
      ggtitle(i)+
      xlim(-1,1)+
      ylim(-1,1)+
      theme(legend.position = "bottom")
  )
}
dev.off()


###### Calculate the percentage of observations in each quarter of combirank scatterplot
# |a|b|
# |c|d|  
coi2 <- c("Aerenchyma", "DMBQ","Syringicacid","Vanillicacid")
coi3 <- c("Aerenchyma", "Suberin")

CountByQuarter <- function(tdf, coi){
  outmat <- matrix(data=NA, ncol=length(coi), nrow=4,
                   dimnames=list(c('A','B','C','D'), coi))
  for(i in coi){
    tlist <- list(a = table((tdf[,i] > 0 & tdf[,"Strigaattachmentspergroot"] < 0)[,1]),
                  b = table((tdf[,i] > 0 & tdf[,"Strigaattachmentspergroot"] > 0)[,1]),
                  c = table((tdf[,i] < 0 & tdf[,"Strigaattachmentspergroot"] < 0)[,1]),
                  d = table((tdf[,i] < 0 & tdf[,"Strigaattachmentspergroot"] > 0)[,1]))
    outmat[,i] <- sapply(tlist, function(x) x['TRUE'])
  }
  return(outmat)
}

quarters2 <- lapply(combo2, CountByQuarter, coi=coi2)
quarters3 <- lapply(combo3, CountByQuarter, coi=coi3)

#If NA, it's 0

##### Calculate contribution of each Phylum 
Phylum <- function(x) {
  x  %>% 
    select(term:Species,contains("combined")) %>% 
    group_by(Phylum) %>%
    summarise_each(funs(sum(!is.na(.)))) %>%
    gather(key = "trait", value = "counts", -(Phylum:Species)) %>%
    group_by(trait) %>%
    mutate(freq = counts/sum(counts))
}

phyl2 <- combo2 %>% lapply (Phylum)
phyl3 <- combo3 %>% lapply (Phylum)


##### Bar plots for the taxa with correlations filtered at 0.2 or -0.2
selcor2_02 <- function(x) {
  x %>% 
    replace_with_na_at(.vars = c("DMBQ", "Strigaattachmentspergroot", "Syringicacid", "Vanillicacid"),
                       condition = ~.x > -0.2) %>%
    replace_with_na_at(.vars = c("Aerenchyma"),
                       condition = ~.x < 0.2)
}

Fun2_list02 <- FunCor2 %>% lapply(selcor2_02)


selcor3_02 <- function(x) {
  x %>% 
    replace_with_na_at(.vars = c("Strigaattachmentspergroot"),
                       condition = ~.x > -0.2) %>%
    replace_with_na_at(.vars = c("Aerenchyma", "Suberin"),
                       condition = ~.x < 0.2)
}

Fun3_list02 <- FunCor3 %>% lapply(selcor3_02)

# Combiranks for taxa with correlations filtered at 0.2 or -0.2 and counts at the phylum level
Fun2_ranks02 <- Fun2_list02 %>% lapply (ranking2)
Fun3_ranks02 <- Fun3_list02 %>% lapply (ranking3)

combo2_02 <- map2(FunCor2, Fun2_ranks02, full_join)
lapply(combo2, head, n=2)

combo3_02 <- map2(FunCor3, Fun3_ranks02, full_join)
lapply(combo3, head, n=2)

phyl2_02 <- combo2_02 %>% lapply (Phylum)
phyl3_02 <- combo3_02 %>% lapply (Phylum)


pdf("Bar plots -combined rank2 - Fungi.pdf")
for(i in names(phyl2_02)){
  print(
    pie <- ggplot(as.data.frame(phyl2_02[[i]]),  aes(x=Phylum, y=counts, width=0.75))+
      geom_bar(width = 1, stat = "identity") +
      ggtitle(i)+
      facet_wrap(~trait, ncol=1)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()


pdf("Bar plots -combined rank3 - Fungi.pdf")
for(i in names(phyl3_02)){
  print(
    pie <- ggplot(as.data.frame(phyl3_02[[i]]),  aes(x=Phylum, y=counts,width=0.75))+
      geom_bar(width = 1, stat = "identity") +
      ggtitle(i)+
      facet_wrap(~trait, ncol=1)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()

Class1 <- function(x) {
  x %>%  
    select(term:Species,contains("combined")) %>% 
    group_by(Phylum, Class) %>%
    summarise_each(funs(sum(!is.na(.)))) %>%
    gather(key = "trait", value = "counts", -(Phylum:Species)) %>%
    filter(str_detect(Class, "^c")) %>%
    filter(Phylum != "punidentified") %>%
    filter(!str_detect(Class, "cunidentified"))
  
}


class2_02 <- combo2_02 %>% lapply (Class1)
class3_02 <- combo3_02 %>% lapply (Class1)

pdf("Bar plots -combined rank3 - Fungi Class .pdf")
for(i in names(class3_02)){
  print(
    pie <- ggplot(as.data.frame(class3_02[[i]]),  aes(x=reorder(Class,desc(Phylum)), y=counts, fill=Phylum, width=0.75))+
      geom_bar(width = 1, stat = "identity", position ="dodge") +
      ggtitle(i)+
      facet_wrap(~trait, ncol=1)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()

pdf("Bar plots -combined rank2 - Fungi Class .pdf")
for(i in names(class2_02)){
  print(
    pie <- ggplot(as.data.frame(class2_02[[i]]),  aes(x=reorder(Class,desc(Phylum)), y=counts, fill=Phylum, width=0.75))+
      geom_bar(width = 1, stat = "identity", position ="dodge") +
      ggtitle(i)+
      facet_wrap(~trait, ncol=1)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()

##### Vann diagrams showing overlap of taxa associated with each traits (only those that passed +/- 0.2 cutoff)

getlistofTermsbyNA <- function(coi, selvec){return(selvec[!is.na(coi)])}
getVennInputList2   <- function(tdf) lapply(tdf[,c("Aer_combined_rank","DMBQ_combined_rank","Syringic_combined_rank","Vanillic_combined_rank")], getlistofTermsbyNA, selvec=tdf$term)
combinedVennFunc2 <- function(df_oi) venn(getVennInputList2(df_oi))
getVennInputList3   <- function(tdf) lapply(tdf[,c("Aer_combined_rank","Sub_combined_rank")], getlistofTermsbyNA, selvec=tdf$term)
combinedVennFunc3 <- function(df_oi) venn(getVennInputList3(df_oi))

pdf("Overlaps-0.2cutoff - Fungi.pdf")
lapply(combo2_02, combinedVennFunc2)
lapply(combo3_02, combinedVennFunc3)
dev.off()

#Extract overlaping taxa btw all traits

# vennObject2 <- lapply(combo2_02, combinedVennFunc2)
vennObject3 <- lapply(combo3_02, combinedVennFunc3)
getIntersect <- function(vennObject){attributes(vennObject)$intersections}

# Intersects2 <- lapply(vennObject2, getIntersect)
Intersects3 <- lapply(vennObject3, getIntersect)

#combine in one df to plot taxonomic membership
AI <- list(
  BK3 = left_join((data.frame(term = Intersects3[["FunBK3week"]][[2]])),combo3[["FunBK3week"]], by="term"),
  Rh3 = left_join((data.frame(term = Intersects3[["FunRh3week"]][[2]])),combo3[["FunRh3week"]], by="term"),
  Endo3 = left_join((data.frame(term = Intersects3[["FunEndo3week"]][[1]])),combo3[["FunEndo3week"]], by="term"),
  EndoSand3 = left_join((data.frame(term = Intersects3[["FunEndoSand3week"]][[1]])),combo3[["FunEndoSand3week"]], by="term"))

AI_phyl <- AI %>% lapply (Phylum)

pdf("Bar plots -Intersects -Fungi.pdf")
for(i in names(AI_phyl)){
  print(
    pie <- ggplot(as.data.frame(AI_phyl[[i]]),  aes(x=Phylum, y=counts, width=0.75))+
      geom_bar(width = 1, stat = "identity") +
      ggtitle(i)+
      facet_wrap(~trait, ncol=1)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()

#Count at the level of Class
Class <- function(x) {
  x   %>%
    select(term:Species,contains("combined")) %>% 
    mutate('name' = paste0(Phylum, '/', Class)) %>%
    count(name) 
}

AI_class <- AI %>% lapply (Class)

pdf("Bar plots -Intersects class - Fungi.pdf")
for(i in names(AI_class)){
  print(
    pie <- ggplot(as.data.frame(AI_class[[i]]),  aes(x=name, y=n, width=0.75))+
      geom_bar(width = 1, stat = "identity") +
      ggtitle(i)+
      scale_y_continuous(position = "right")+
      theme_minimal()+
      theme(
        axis.text.x=element_text(angle = -90, hjust = 0)))
  
  
}  
dev.off()

#Barplot for overlap of suberin and aerenchyma with all compartments
bk3 <- AI[["BK3"]] %>% select(term:Species) %>% mutate(comp = "Bulk soil") %>% mutate('name' = paste0(Phylum, '/', Class))
rh3 <- AI[["Rh3"]] %>% select(term:Species) %>% mutate(comp = "Rhizosphere") %>% mutate('name' = paste0(Phylum, '/', Class))
en3 <- AI[["Endo3"]] %>% select(term:Species) %>% mutate(comp = "Soil-associated endosphere") %>% mutate('name' = paste0(Phylum, '/', Class))
ensand3 <- AI[["EndoSand3"]] %>% select(term:Species) %>% mutate(comp = "Sands-associated endosphere") %>% mutate('name' = paste0(Phylum, '/', Class))

wpi3 <- full_join((full_join(bk3,rh3)),(full_join(en3,ensand3))) %>% count(comp,name) 

pdf("Bar plots -Intersects wpi3 class fungi.pdf")
ggplot(wpi3,  aes(x=name, y=n, fill=comp,  width=0.75))+
  geom_bar(width = 1, stat = "identity", position ="dodge") +
  scale_y_continuous(position = "right")+
  theme_minimal()+
  theme(
    axis.text.x=element_text(angle = -90, hjust = 0))
dev.off()


###Create tables with all ranks per trait

tof <- c("Aer_combined_rank", "DMBQ_combined_rank", "Syringic_combined_rank", "Vanillic_combined_rank", "Sub_combined_rank")
contains_tof <- function(x, nam) sum(grepl(nam, colnames(x))) == 1
sel_names <- c("term", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tofL <- vector(length = length(tof), mode='list')
names(tofL) <- tof
for(i in tof){
  Fun <- Fun_ranks[sapply(Fun_ranks, contains_tof, nam = i)]
  t.sel_names <- c(sel_names, i)
  trait <- lapply(Fun, "[", t.sel_names)
  
  trait2 <- lapply(1:length(trait), function(x) {
    tL <- trait [[x]]
    colnames(tL) <- c("term","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", names(trait)[x])
    return(tL)})
  
  tofL[[i]] <- Reduce(
    function(x, y, ...) full_join(x, y, by = c('term',"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
    trait2)
}
lapply(tofL, head, n=2)

write.xlsx(tofL, "Fungi Combined ranks trait*compartment.xlsx")


