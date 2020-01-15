#!/usr/bin/env Rscript

args <- commandArgs(trailing=T)

if(len(args)==0) {
    cat("\n  usage: makeMPG.R origMappingFile pooledNormalMappingFile manifestFile\n\n")
    quit()
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))

origMappingFile=args[1]
pooledNormalMappingFile=args[2]

fixSampleName<-function(x){
    basename(x) %>%
        gsub("Sample_","s_",.) %>%
        gsub("_IGO_.*","",.) %>%
        gsub("-","_",.)
}

mapping=read_tsv(origMappingFile,col_names=F) %>%
    mutate(X2=fixSampleName(X4)) %>%
    arrange(factor(X2,levels=mixedsort(unique(X2)))) %>%
    mutate(IID=gsub(".*Sample_","",X4)) %>%
    mutate(IID=gsub("-","_",IID)) %>%
    select(-X4,X4)

manifest=read_xlsx(args[3]) %>%
    mutate_all(~gsub("-","_",.)) %>%
    mutate_all(~gsub("[^A-Za-z0-9_]","",.))
colnames(manifest)[1]="FID"

numDuplicateSamples=manifest %>% count(FID) %>% filter(n>1) %>% nrow
if(numDuplicateSamples>0) {
    cat("\n\nError in manifest; duplicate samples\n")
    ss=knitr::kable(manifest %>% count(FID) %>% filter(n>1))
    cat(paste(ss,collapse="\n"))
    cat("\n\n")
    stop("ERROR")
}

mapping=left_join(mapping,manifest,by=c(IID="FID"))

#
# Make a neuron pool for each patient
#

patients <- mapping %>% distinct(Patient) %>% pull

map.pool=map(
        patients,
        ~ filter(mapping,Patient == . & Type=="NORMAL") %>%
            mutate(X2=cc("s",Patient,"NEUN","POOL")) %>%
            mutate(Type="NEUN.POOL")
        ) %>%
    bind_rows


mapping=bind_rows(mapping,map.pool)

map.normpool=read_tsv(pooledNormalMappingFile,col_names=F)

if(nrow(map.normpool)>0) {

    map.normpool = map.normpool %>%
        mutate(X2=fixSampleName(X4)) %>%
        mutate(Type="POOLED.NORMAL") %>%
        mutate(Patient="CTRL.POOL")

    mapping=bind_rows(mapping,map.normpool)

}

#
# Pair to NUEN.POOL
#

getNEUNPoolPair<-function(pii) {
    neunPool=mapping %>% filter(Patient==pii & Type=="NEUN.POOL") %>% distinct(X2) %>% pull
    if(len(neunPool)>0) {
        mapping %>%
            filter(Patient==pii & Type=="TUMOR")  %>%
            select(X2) %>%
            mutate(Normal=neunPool) %>%
            select(Normal,X2)
    } else {
        return(NULL)
    }
}

pair.neun=map(patients,getNEUNPoolPair) %>% bind_rows

if(nrow(map.normpool)>0) {
    pair.pool=mapping %>%
        filter(!grepl("POOL",Type)) %>%
        distinct(X2) %>%
        mutate(Normal="s_FROZENPOOLEDNORMAL") %>%
        select(Normal,X2)

    pair=bind_rows(pair.neun,pair.pool) %>% distinct
} else {
    pair = distinct(pair.neun)
}

projectTag=basename(getwd())

write_tsv(mapping %>% select(X1,X2,X3,X4,X5),cc(projectTag,"sample_mapping.txt"),col_names=F)
write_tsv(pair,cc(projectTag,"sample_pairing.txt"),col_names=F)

grouping <- mapping %>%
    select(X2,Patient) %>%
    mutate(GID=cc("Group",as.numeric(factor(Patient)))) %>%
    select(-Patient) %>%
    arrange(GID,X2) %>%
    distinct

gCounts=count(grouping,GID) %>% pull(n)
#
# Force new groupoing convention
#
if(any(gCounts>0)) {
    grouping=mapping %>%
        select(X2) %>%
        distinct %>%
        mutate(Group=cc("Group",row_number()))
}

write_tsv(grouping,cc(projectTag,"sample_grouping.txt"),col_names=F)
