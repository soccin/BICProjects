#!/usr/bin/env Rscript

args <- commandArgs(trailing=T)

if(len(args)==0) {
    cat("\n  usage: findMousePools.R mappingFile\n\n")
    quit()
}

require(readr)
require(dplyr)
require(purrr)
require(fs)
args=commandArgs(trailing=T)

mapping=read_tsv(args[1],col_names=F)

paths=mapping %>%
    mutate(SeqDir=gsub(".Project.*","",X4)) %>%
    distinct(SeqDir) %>%
    map(dir_ls,regexp="Project_POOLED") %>%
    map(dir_ls,recurse=T,regexp="MOUSEPOOLEDNORMAL") %>%
    unlist %>%
    as.character

#write(paths[dir.exists(paths)],"sampleDirs.pool")

newmap=tibble(X4=paths[dir.exists(paths)]) %>%
    mutate(X1="_1",X5="PE") %>%
    mutate(X2=gsub(".*Sample_","s_",X4) %>% gsub("_IGO_.*","",.)) %>%
    mutate(X3=gsub(".Project_.*","",X4) %>% gsub(".*FASTQ.","",.)) %>%
    select(X1,X2,X3,X4,X5)

write_tsv(newmap,"mapping.pool",col_names=F)
