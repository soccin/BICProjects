#!/ifs/work/socci/opt/R/3.4.3/bin/Rscript --no-save

args <- commandArgs(trailing=T)

if(len(args)==0) {
    cat("\n  usage: fixMappingIDs.R mappingFile\n\n")
    quit()
}

require(readr)
require(dplyr)
args=commandArgs(trailing=T)

map=read_tsv(args[1],col_names=F)
map1=map %>% mutate(X2=cc("s",gsub(".*Sample_","",X4) %>% gsub("_IGO_.*","",.)))

write_tsv(map1,"mapping.fix",col_names=F)
