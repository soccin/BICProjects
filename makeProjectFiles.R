#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(yaml))
suppressPackageStartupMessages(require(stringr))
cArgs=commandArgs(trailing=T)

len <- function(x) {length(x)}
cc <- function(...) {paste(...,sep='_')}


#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#

# Set defaults first

cat("\n\nmakeProjectFiles.R Version 1.0.1\n\n")

args=list(ASSAY=NULL)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(purrr))
suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(yaml))

pwd=getwd()

projNo=gsub("Proj_","",basename(pwd))
workflow=basename(dirname(pwd))

readmeFile=dir_ls(regexp="README")
if(len(readmeFile)>0) {
    readmeDat=readLines(readmeFile) %>% gsub("^ *","",.)

    if(len(grep(" - ",readmeDat))>3) {

        # Old/Orig style of README from #REQUEST forms
        readme=strsplit(readmeDat," - ") %>% map(2)
        names(readme)=strsplit(readmeDat," - ") %>% map(1) %>% unlist %>% make.names

    } else {

        # New style request

        readme=read_yaml(readmeFile)
        names(readme)=gsub("[ ?/]",".",names(readme))

    }

} else {
    cat("\n   No README File. Need some basic info (Cost center/Fund number)\n\n")
    quit()
}

drafts=file.path("/juno/projects/BIC/drafts",cc("Proj",projNo))

#
# If no mapping file get the one from drafts
#
mappingFile=dir_ls('.',regexp="_sample_mapping.txt$")
if(len(mappingFile)!=1) {
    xx=file.copy(dir_ls(drafts,regexp="_sample_mapping.txt"),pwd)
    if(len(xx)<1 || !xx) {
        cat("\n    ERROR: Can not file mapping file\n\n")
        quit()
    }
}

if(file.exists(drafts) && len(dir_ls(drafts,regexp="_request.txt"))>0) {
    requestDat=readLines(dir_ls(drafts,regexp="_request.txt"))
    request=strsplit(requestDat,": ") %>% map(2)
    names(request)=strsplit(requestDat,": ") %>% map(1) %>% unlist
} else {
    request=list()
}

if(file.exists("_request")) {
    requestDat=readLines("_request")
    request0=strsplit(requestDat,": ") %>% map(2)
    names(request0)=strsplit(requestDat,": ") %>% map(1) %>% unlist
    request[names(request0)]=request0
}

if(file.exists("_request.txt")) {
    requestDat=readLines("_request.txt")
    request0=strsplit(requestDat,": ") %>% map(2)
    names(request0)=strsplit(requestDat,": ") %>% map(1) %>% unlist
    request[names(request0)]=request0
}

# Turn this off 2021-08-16
# #######################################################################################
# # Use investigator info from README
# #######################################################################################
# if(!is.null(readme$Your.name)) {
#     request$Investigator_Name=stringr::str_to_title(readme$Your.name)
# }

# if(!is.null(readme$Your.email)) {
#     request$`Investigator_E-mail`=readme$Your.email
#     request$Investigator=gsub("@.*","",readme$Your.email)
# }

#######################################################################################
# Workflow specific variables
#######################################################################################

rnaSeqDifferential=FALSE
if(workflow == "rnaseq") {

    cat("RNAseq workflow processing ...")

    if(readme$What.analysis.are.you.requesting.=="Full analysis/differential gene lists") {

        rnaSeqDifferential=TRUE

    }

    if(request$LIMS_Strand=="stranded-reverse") {
        request$Strand="Reverse"
    } else if(request$LIMS_Strand=="non-stranded") {
        request$Strand="None"
    } else {
        cat("\n\n  UNKNOWN STRAND VALUE CHECK MANUALLY\n\n")
        quit("ERROR::UNKNOWN STRAND")
    }

    if(rnaSeqDifferential) {
        request$`Charges-Service`="RNASeq_Gene_Differential"
    } else {
        request$`Charges-Service`="RNASeq_Gene_Counts"
    }



    cat("\n")

}

request$Pipelines=case_when(
    workflow == "variant" ~ "variants",
    workflow == "chipseq" ~ "ChIP-seq Mapping",
    workflow == "other" ~ "ChIP-seq Mapping",
    workflow == "rnaseq" ~ ifelse(rnaSeqDifferential,
                                    "RNASEQ_DIFFERENTIAL_GENE_V1",
                                    "RNASEQ_STANDARD_GENE_V1"
                                  ),
    T ~ "NA"
    )

request$Run_Pipeline=case_when(
    workflow == "variant" ~ "variants",
    workflow == "chipseq" ~ "chipseq",
    workflow == "other" ~ "chipseq",
    workflow == "rnaseq" ~ "rnaseq",
    T ~ "NA"
    )

if(is.null(request$Assay)) {
    request$Assay=case_when(
        workflow == "chipseq" ~ "na",
        workflow == "other" ~ "na",
        T ~ "na"
        )
}

request=request[request!="NA"]
request=request[!(map(request,is.null) %>% unlist)]
request$ProjectFolder=pwd

#
# Try to fill out known values from README
#

assayInfoPath="/ifs/projects/BIC/bin/data"
targetAssayFld=grep("Target.assay",names(readme),value=T)
if(len(targetAssayFld)>0) {
    targetAssayName=readme[[targetAssayFld]]
    assayYaml=file.path(assayInfoPath,paste0(targetAssayName,".yaml"))
    if(file.exists(assayYaml)) {
        assayInfo=read_yaml(assayYaml)
        for(ni in names(assayInfo)) {
            request[[ni]]=assayInfo[[ni]]
        }
    }
}

#
# Get necessary values if they have not already been set
#

if(is.null(request$NumberOfSamples)) {
    suppressPackageStartupMessages(library(readr))
    request$NumberOfSamples <- read_tsv(mappingFile,col_names=F) %>% distinct(X2) %>% nrow
}

if(is.null(request$ProjectID)) {
    request$ProjectID <- basename(request$ProjectFolder)
}

if(is.null(request$RunNumber)) {
    request$RunNumber <- 1
}

if(is.null(request$`PI_E-mail`)) {
    request$`PI_E-mail` <- readme$Lab.head...PI.email
    request$PI <- gsub("@.*","",request$`PI_E-mail`)
}

if(is.null(request$`Investigator_E-mail`)) {
    request$`Investigator_E-mail` <- readme$Your.email
    request$Investigator <- gsub("@.*","",request$`Investigator_E-mail`)
}

if(is.null(request$PI)) {
    request$PI <- gsub("@.*","",request$`PI_E-mail`)
}

if(is.null(request$Investigator)) {
    request$Investigator <- gsub("@.*","",request$`Investigator_E-mail`)
}

if(is.null(request$DeliverTo_Name)) {
    request$DeliverTo_Name=stringr::str_to_title(readme$Your.name)
    request$DeliverTo_Email=readme$Your.email
}

# if(is.null(request$PI_Name)) {
#     request$PI_Name=""
# }

# if(is.null(request$Investigator_Name)) {
#     request$Investigator_Name=""
# }

if(is.null(request$Species)) {
    request$Species <- readme$Species
}

# if(is.null(request$)) {
# }

#
# Validate
#

if(workflow=="variant") {
    if(is.null(request$Assay)) {
        cat("\n\tCan not have NULL Assay field in variant workflow\n\n")
        quit()
    }
    if(is.null(request$Comments)) {
        request$Comments="Run additional variant callers - muTect2, Vardict, Strelka"
    }
}

newRequestFile=cc("Proj",projNo,"request.txt")

chargeFlds=grep("Charges-",names(request))
if(len(chargeFlds)>0) {
    request1=request[-chargeFlds]
} else {
    request1=request
}

write(paste0(names(request1),": ",unlist(request1)),newRequestFile)
cat("\n",file=newRequestFile,append=TRUE)

if(is.null(request[["Charges-Division"]])) {
    request[["Charges-Division"]]="BIC"
}

ccfn=paste0(readme$Cost.center,"/",readme$Fund.number)
cat("Charges-CCFN:",ccfn,"\n",file=newRequestFile,append=TRUE)
cat("Charges-Division:",request[["Charges-Division"]],"\n",file=newRequestFile,append=TRUE)
cat("Charges-ProjectNumber:",projNo,"\n",file=newRequestFile,append=TRUE)
cat("Charges-Qty:",request$NumberOfSamples,"\n",file=newRequestFile,append=TRUE)
cat("Charges-Service:",request[["Charges-Service"]],"\n\n",file=newRequestFile,append=TRUE)

# if(is.null(request[["Charges-Service"]])) {
#     cat("\n   Need to added Charges-Service code in [_request] file\n\n")
#     cat("Charges-Service: \n",file="_request",append=TRUE)
# }

