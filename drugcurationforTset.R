getwd()
library(dplyr)
install.packages("gsubfn")
library(gsubfn)
install.packages("textclean")
library(textclean)

dixa <- read.csv(file = "s_Hepatocyte.csv", header = TRUE, sep = "\t")
lab <- read.csv(file = "drugs_with_ids.csv", header = TRUE)


labdruguniq <- lab$unique.drugid

as.data.frame(labdruguniq)

dixadrug <- dixa$Factor.Value.Compound.
ncol(dixadrug)
dixadrug

dixadruguniq <- as.data.frame(unique(dixadrug))

#below dataframes contain list of drugs. dixadrug has 124 relevant compounds, labdruguniq has thousands, we have to 
#isolate the relevant compounds and map to dixadruguniq 
dixadruguniq
labdruguniq

intersect(labdruguniq, dixadruguniq)
#no intersect

dixadruguniq <- mutate_all(dixadruguniq, list(toupper))

dixadruguniq[,order(colnames(dixadruguniq))]
#alphabetically ordered 
dixadruguniq

dixadruguniq$`unique(dixadrug)` <- gsub('\\s+', '', dixadruguniq$`unique(dixadrug)`)
dixadruguniq



dixadruguniqA <- as.data.frame(mgsub(dixadruguniq$`unique(dixadrug)`, c("-", "'", ","), c("")))
dixadruguniqA


#------------

labdruguniq <- as.data.frame(unique(labdruguniq))
labdruguniq <- mutate_all(labdruguniq, list(toupper))
is.data.frame(labdruguniq)

labdruguniq

labdruguniq[,order(colnames(labdruguniq))]


labdruguniqA <- as.data.frame(mgsub(labdruguniq$`unique(labdruguniq)`, c("-", "'", ",", "[", "]", ":", ";", "(", ")", ""), c("")))
labdruguniqA

head(labdruguniqA)

names(dixadruguniqA) <- c("drug_id")
names(labdruguniqA) <- c("drug_id")

labdixacommon <- intersect(labdruguniqA,dixadruguniqA)

dixadruguniqA

labdixacommon
difference <- setdiff(dixadruguniqA, labdixacommon)
        
difference[,order(colnames(difference))]

difference
is.data.frame(difference)        

labdixacommon

#drugcuration. two columns. one for unique identifier (lab annotated drug names), second with column from the drugs in the study 


drug_id <- c('ANTIMYCIN A', 'ESTRADIOL', 'CADMIUM DICHLORIDE', 'N7S12310TETRAMETHOXY9OXO5679TETRAHYDROBENZOAHEPTALEN7YLACETAMIDE', 'CYCLOSPORIN A', 'CYTOCHALASIN B', 'ZONALON', 'GENTAMICIN','KANAMYCIN','XOLEGEL','LITHOCHOLIC ACID','TANDUTINIB')

missedcommon <- data.frame(drug_id)

total <- rbind(missedcommon,labdixacommon)
total

labdruguniqA
labdruguniqA[labdruguniqA$drug_id == 'ANTIMYCIN A',]

nrow(labdruguniqA)
labdruguniqA$drug_id[48755]
labdruguniqA$drug_id[189]
labdruguniqA$drug_id[48203] == 'DOXEPIN'
labdruguniqA$drug_id[49348]

test <- intersect(labdruguniqA, total)
test

labdixacommon
setdiff(total,test)

total

