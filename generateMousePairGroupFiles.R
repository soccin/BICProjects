require(readr)
require(dplyr)
require(fs)
require(stringdist)

mapFile=dir_ls(regex="_sample_mapping.txt")
map=read_tsv(mapFile,col_names=F)
samples=map %>% distinct(X2) %>% arrange(X2) %>% pull

manifest=read_tsv(dir_ls(regex="_metadata_samples.txt"))

normalIds=manifest %>%
    filter(tumorOrNormal=="Normal") %>%
    pull(investigatorSampleId) %>%
    gsub("-","_",.) %>%
    cc("s",.)

poolIds=grep("POOLEDNORMAL",samples,value=T)

tumorIds=setdiff(samples,c(normalIds,poolIds))

pairs=list()
for(ni in normalIds) {
    tii=which.min(stringdist(ni,tumorIds))
    pairs[[len(pairs)+1]]=tibble(Normal=ni,Tumor=tumorIds[tii])
    pairs[[len(pairs)+1]]=tibble(Normal=poolIds,Tumor=ni)
}
pairs=bind_rows(pairs) %>% arrange(Normal)

write_tsv(pairs,gsub("_sample_mapping.txt","_sample_pairing.txt",mapFile),col_names=F)

map %>%
    distinct(X2) %>%
    mutate(Group=cc("Group",row_number())) %>%
    write_tsv(gsub("_sample_mapping.txt","_sample_grouping.txt",mapFile),col_names=F)
