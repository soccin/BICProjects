#!/usr/bin/env Rscript

args <- commandArgs(trailing=T)

if(len(args)==0) {
    cat("\n  usage: fixMappingIDs.R mappingFile\n\n")
    quit()
}

suppressPackageStartupMessages(require(readr))
suppressPackageStartupMessages(require(dplyr))
args=commandArgs(trailing=T)

map=read_tsv(args[1],col_names=F)
map1=map %>% mutate(X2=cc("s",gsub(".*Sample_","",X4) %>% gsub("_IGO_.*","",.)))

write_tsv(map1,"mapping.fix",col_names=F)
