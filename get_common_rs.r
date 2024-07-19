# Load necessary libraries
library(dplyr)
library(tableone)
library(data.table)
library(tidyverse)


# ADNI3 genotyping data available?
t1 = fread('PLINK_set1/ADNI3_PLINK_Final.bim')
t2 = fread('PLINK_set2/ADNI3_PLINK_FINAL_2nd.bim')
# t1 %>% filter(V1==8) %>% head(20)
t1s = t1 %>% filter(V5 %in% c('A', 'T', 'C', 'G')) %>% filter(V6 %in% c('A', 'T', 'C', 'G')) %>% mutate(newID = paste0('chr', V1, '-', V4, ':', V2, '-', V5, '..', V6))
t2s = t2 %>% filter(V5 %in% c('A', 'T', 'C', 'G')) %>% filter(V6 %in% c('A', 'T', 'C', 'G')) %>% mutate(newID = paste0('chr', V1, '-', V4, ':', V2, '-', V5, '..', V6))
commonVar = intersect(t1s$newID, t2s$newID) # a,t,g,c are also considered
# write out the list of commonRS
commonVar2 = t1s %>% filter(newID %in% commonVar) %>% distinct(V2, .keep_all=T) %>% distinct(V1, V4, .keep_all=T) %>% select(V2)
commonVar2 %>% fwrite('adni3_common_var.txt', col.names=F)
commonVar2 = commonVar2$V2
print('Number of common Variants')
print(length(commonVar2))


# RSID table
commonRS = unique(grep('rs', commonVar2, value=T))
dt = data.frame(ID = commonRS) %>% 
    mutate(ID2 = gsub('^GSA-', '', ID)) %>%
    mutate(ID2 = gsub('^exm-', '', ID2)) %>%
    mutate(ID2 = gsub('^new', '', ID2)) %>%
    separate(ID2, c('ID3', 'ID4'), sep=':', fill='left', remove=F) %>% 
    separate(ID4, c('ID5', 'ID6'), sep='-', fill='right', remove=F) %>% 
    rename(rsid = ID5) %>%
    filter(ID!=rsid) %>%
    select(ID, rsid)

fwrite(dt, 'adni3_common_ID_rsID_mapping.txt', col.names=F, quote=F, sep='\t')
print('ID<-->RSID mapping table written')
dim(dt)
print('Done')