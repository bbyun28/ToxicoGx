
phenodataa <- read.csv(file = "/Desktop/hepatocytemicroarray/s_Hepatocyte.csv", sep = "\t")

library(ToxicoGx)
library(Biobase)
library(dplyr)
#___________________

##PHENODATA##

phenodata <- read.csv(file = "s_Hepatocyte.csv", head = TRUE, sep="\t")
phenodata
is.data.frame(phenodata)

#cleaning up phenoData table 
phenodata[phenodata == 'NA' | phenodata =='not specified'] <- NA
phenodataC <- phenodata[,colSums(is.na(phenodata)) < nrow(phenodata)]
ncol(phenodata)
ncol(phenodataC)
names(phenodataC) <- gsub("Characteristics.", "", names(phenodataC))
phenodataC$DoseUnit. <- "μM"
names(phenodataC) <- gsub("Factor.Value", "", names(phenodataC))
phenodataC$Organism. <- "R.norvegicus"
names(phenodataC) <- gsub("Term.", "", names(phenodataC))
phenodataC <- phenodataC[,-32]
ncol(phenodataC)
colnames(phenodataC)[colnames(phenodataC) == "CompoDosetDoseReplicate#"] <- "Compound.DoseDuration.Dose.BioReplicate#"
names(phenodataC) <- gsub("Strain.", "Strain", names(phenodataC))
names(phenodataC) <- gsub("Organism.", "Organism", names(phenodataC))
names(phenodataC) <- gsub("Subject.ID.", "Subject.ID", names(phenodataC))
names(phenodataC) <- gsub("Sex.", "Sex", names(phenodataC))
names(phenodataC) <- gsub("Cell.", "Cell", names(phenodataC))
names(phenodataC) <- gsub("Assay.Type.", "Assay.Type", names(phenodataC))
names(phenodataC) <- gsub("Biological.Replicate.", "Biological.Replicate", names(phenodataC))
names(phenodataC) <- gsub("Technical.Replicate.", "Technical.Replicate", names(phenodataC))
names(phenodataC) <- gsub(".Compound.", "Compound.Name", names(phenodataC))
names(phenodataC) <- gsub("Control.", "Control", names(phenodataC))
names(phenodataC) <- gsub("Sample.Match.", "Sample.Match", names(phenodataC))
names(phenodataC) <- gsub(".Dose.", "Dose", names(phenodataC))
names(phenodataC) <- gsub("StdInChIKey.", "StdInChIKey", names(phenodataC))
names(phenodataC) <- gsub("Comment.chEMBL.ID.", "Comment.chEMBL.ID", names(phenodataC))
names(phenodataC) <- gsub("DoseUnit.", "Dose.Unit", names(phenodataC))
names(phenodataC) <- gsub("DoseDuration.", "Dose.Duration", names(phenodataC))
names(phenodataC) <- gsub("Dose.DurationUnit.", "Dose.Duration.Unit", names(phenodataC))
names(phenodataC) <- gsub("Vehicle.", "Vehicle", names(phenodataC))
colnames(phenodataC)[colnames(phenodataC)=="Source.REF"] <- "Org.ID.Abbr"
phenodataC$Org.ID.Abbr <- "NCBIT"
phenodataC$Source.REF.5 <- "CHEBI"
colnames(phenodataC)[colnames(phenodataC)=="Source.REF.5"] <- "CompoundName.abbr"
colnames(phenodataC)[colnames(phenodataC)=="Accession.Number.5"] <- "Acccesion.Number"
phenodataC$Accession.Number <- phenodataC$Protocol.REF <- phenodataC$Source.REF.4 <- phenodataC$Accession.Number.4 <- phenodataC$Source.REF.6 <- NULL
phenodataC$Source.REF.7 <- phenodataC$Protocol.REF <- phenodataC$Accession.Number.6 <- phenodataC$Accession.Number <- phenodataC$StdInChIKey <- phenodataC$Acccesion.Number <- NULL
colnames(phenodataC)[colnames(phenodataC)=="Dose.Duration"] <- "Dose.Duration.inDays"


ncol(phenodataC)
nrow(phenodataC)

ncol(phenodataC)
phenodataC[827:939,] 
getOption("max.print", default=NULL)

phenodataC$Dose.Duration.inDays
unique(phenodataC$Dose.Duration.inDays)
#112 controls in total, 52 and 61 
#0.67, 1.00 in terms of dose duration in days 

#check how many times each drug repeats, then we know how many replicates we have
#check within those replicates, how many times each is in 0.67 days and 1.0 days, respectively
#see if 113 controls match (should have one control in each set of 0.67 and 1.0 days replicates)


phenodataC$Sample.Match

unique(phenodataC$Compound.Name)
#125 drugs 

phenodataC

phenodataC$Source.Name[827:939]
Ctrlonly <- phenodataC$Subject.ID[827:939]
Rest <- phenodataC$Subject.ID[1:826]

RestDose <- phenodataC$Dose.Duration.inDays[1:826]
unique(RestDose)
sort(RestDose, decreasing = FALSE,)



Ctrlonly
intersect(Ctrlonly, Rest)


RestDoseData <-as.data.frame(RestDose)
is.data.frame(RestDoseData)
count(RestDoseData, vars = "1.00")

table(phenodataC$Compound.Name)


phenodataC$Source.Name
#113 controls in total 


freq <-ave(rep(1, times=nrow(phenodataC)), phenodataC$Compound.Name, FUN=sum)
            phenodataC[order(freq,phenodataC$Compound.Name),]


saveRDS(phenodataC, file = "/Users/parwaiznijrabi/Desktop/phenoData.rds")

