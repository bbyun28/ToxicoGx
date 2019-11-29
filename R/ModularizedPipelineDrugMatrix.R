library(ToxicoGx)
library(Biobase)
library(affy)
library(affyio)
library(BiocManager)
library(rat2302rnensgcdf)
library(dplyr)
library(biomaRt)

#install.packages("affy")
#install.packages("affyio")
install.packages("BiocManager")
#install.packages("rat2302rnensgcdf")
#install.packages("~/Desktop/rat2302rnensgcdf", repos= NULL, type="source")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")


devtools::install_github("bhklab/ToxicoGx", ref = "master")
library(ToxicoGx)

downloadTSet('drugMatrix')


####PHENODATA####
#s_hepatocyte file can be found on drugMatrix dixa database. metadata for hepatocyte data 
#939 relevant CEL files from 5590 total file set - the microarray data for all tissue types from drugMatrix 
#cleaning up phenoData table

#changes to original phenodata file 
#colname changes, multiple colname changes
#removing all NA columns
#adding chiptype, xptype, organid, batchid, celfilename, duration and other columns 

create_phenoData <- function(species=c("Rat"), verbose = TRUE) { 
  if (verbose) {message("Creating phenoData object...")} 
  s_Hepatocyte <- read.csv("./data/phenodataDMSOcontrols.csv", header = TRUE, stringsAsFactors = FALSE)
  rownames(s_Hepatocyte) <- s_Hepatocyte$samplename
  if(verbose) {message("phenoData object created!")}
  return(s_Hepatocyte)
}

create_phenoData("Rat")

create_exprsData <- function(species=c("Rat"), phenoData, verbose = TRUE) {
  if (verbose) {message("Creating eset object")}

  if(species == "Rat"){
    install.packages("rat2302rnensgcdf")
    library("rat2302rnensgcdf")
    cdf <- "rat2302rnensgcdf"
  }
  esetNorm <- just.rma(filenames = celFiles, verbose = TRUE, cdfname = "rat2302rnensgcdf")
  eset <- esetNorm
  storageMode(eset)<-"environment"
  eset <-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  storageMode(eset)<-"lockedEnvironment"
  annotation(eset)<-"rna"
  if (verbose) {message("eset object created!")}
  return(eset)
}
#celFiles <- list.celfiles("/Users/parwaiznijrabi/Desktop/Hep939", full.names = TRUE)
#eset <- readRDS("./esetNorm.rds")
#hep939 folder contained all cel files relevant to rat hepatocyte. Can be found by checking phenodata
#saveRDS(esetNorm, file = "/Users/parwaiznijrabi/Desktop/esetNorm.rds")  


create_featureData <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating featureData object...")}
  if (species == "Rat"){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)
    storageMode(eset) <- "environment"
    affxrows <- rownames(eset@assayData$exprs)
    rownames(eset@assayData$exprs) <- substr(rownames(eset@assayData$exprs), 1, nchar(affxrows)-3)

    saveRDS(affxrows, file = "./data/probesRatVitro1.rds") 
    
    CELgenes <- readRDS("./data/probesRatVitro1.rds")
    CELgenes1 <- gsub(".at", " ", CELgenes)
    results <-getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes1, mart=ensembl,checkFilters = TRUE)
    uniqueB <- results[!duplicated(results$ensembl_gene_id),]
    CELnotB <- unique(CELgenes1) [!unique(CELgenes1) %in% uniqueB$ensembl_gene_id]
    names(uniqueB) <- c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
    finalFeature <- uniqueB
    
    finalFeature$BEST <- NA
    names(finalFeature) <- c("Symbol", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id", "BEST")
    rownames(finalFeature) <- finalFeature$gene_id
    finalFeature$gene_id
    geneid1 <- finalFeature$gene_id
    
    for (i in 1:length(geneid1)) {
      geneid1[i] = paste(geneid1[i], "at", sep="_")
    }
    geneid1
    finalFeature$gene_id <- geneid1
    finalFeature$gene_id
    finalFeature[,1]
    rownames(finalFeature) = finalFeature$gene_id
    
    if(verbose) {message("featureData object created!")}
    return(finalFeature)
    
  }
}

#LABANNOT <- read.csv("/Users/parwaiznijrabi/Desktop/ToxicoGx3/TGX/metadata/annot_ensembl_all_genes.csv", header = T)
#CELnotB <- unique(CELgenes1) [!unique(CELgenes1) %in% uniqueB$ensembl_gene_id]
#newLabAnnot <- subset(LABANNOT, LABANNOT$X %in% CELnotB, select=c(X,gene_biotype, gene_id, gene_name, transcript_id, transcript_name,EntrezGene.ID))
#putting together eset


create_Expressionset <- function(species=c("Rat"), eset, verbose = TRUE){
  if (verbose) {message("Creating expressionset...")}
  if (species == "Rat"){
phenoData <- readRDS("./phenodataDMSOcontrols.rds")
featureData <- readRDS("./featureData.rds")
pData(eset)<-phenoData
fData(eset)<-featureData
storageMode(eset)<-"lockedEnvironment"
saveRDS(eset1, file="/Users/parwaiznijrabi/Desktop/ExpressionSetFINALerror.rds")

  }
}


######################toxicoset constructor function ######################

tSet <- ToxicoSet("drugmatrix_hepatocyte",
                  molecularProfiles=list("rna"=ExpressionSet),
                  cell=cell,
                  drug=drug,
                  sensitivityInfo=NULL,
                  sensitivityRaw=NULL,
                  sensitivityProfiles=NULL,
                  curationDrug=curationDrug,
                  curationCell=curationCell,
                  curationTissue=curationTissue,
                  datasetType=c("perturbation"),
                  verify = TRUE)


#toxicosetconstructor 

#all(rownames(fData(esetFINAL))%in% rownames(exprs(esetFINAL))) 
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#length(intersect(pData(esetFINAL)[,"cellid"] , cell$cellid))
#View(esetFINAL@phenoData@data)
#View(esetFINAL@assayData$exprs)
#View(esetFINAL@featureData@data)

# TESTS #
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(pData(eset)[,"cellid"] , cell$cellid))
#length(intersect(unique(sensitivityInfo$cellid) , cell$cellid))