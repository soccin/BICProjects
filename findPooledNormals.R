#!/usr/bin/env Rscript

args <- commandArgs(trailing=T)

if(len(args)==0) {
    cat("\n  usage: findMousePools.R mappingFile [POOL_TYPE]\n\n")
    quit()
}

suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(purrr))
suppressPackageStartupMessages(require(fs))
args=commandArgs(trailing=T)

mapping=read_tsv(args[1],col_names=F,col_types = cols(.default = "c"))

paths01=mapping %>%
    mutate(SeqDir=gsub(".Project.*","",X4)) %>%
    distinct(SeqDir) %>%
    map(dir_ls,regexp="Project_POOLED") %>%
    map(dir_ls)

if(len(args)==1) {

    poolTypes=paths01 %>%
        map(~gsub(".*Sample_","",.)) %>%
        map(~gsub("_[^_]+$","",.)) %>%
        unlist %>%
        as.character %>%
        unique %>%
        sort
    cat("\nPool Types:\n")
    for(pp in poolTypes) {
        cat("  ",pp,"\n")
    }

    cat("\n")

    quit()

} else {

    paths=grep(args[2],unlist(paths01)%>%as.character,value=T)

}

newmap=tibble(X4=paths[dir.exists(paths)]) %>%
    mutate(X1="_1",X5="PE") %>%
    mutate(X2=gsub(".*Sample_","s_",X4) %>% gsub("_IGO_.*","",.)) %>%
    mutate(X3=gsub(".Project_.*","",X4) %>% gsub(".*FASTQ.","",.)) %>%
    select(X1,X2,X3,X4,X5)

write_tsv(newmap,"mapping.pool",col_names=F)
