% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ToxicoSet.R
\docType{class}
\name{ToxicoSet-class}
\alias{ToxicoSet-class}
\alias{.ToxicoSet}
\alias{cellInfo,ToxicoSet-method}
\alias{cellInfo<-,ToxicoSet,data.frame-method}
\alias{drugInfo,ToxicoSet-method}
\alias{drugInfo<-,ToxicoSet,data.frame-method}
\alias{phenoInfo,ToxicoSet,character-method}
\alias{phenoInfo<-,ToxicoSet,character,data.frame-method}
\alias{molecularProfiles,ToxicoSet-method}
\alias{molecularProfiles<-,ToxicoSet,character,character,matrix-method}
\alias{molecularProfiles<-,ToxicoSet,character,missing,matrix-method}
\alias{molecularProfilesSlot,ToxicoSet-method}
\alias{molecularProfilesSlot<-,ToxicoSet,ANY-method}
\alias{featureInfo,ToxicoSet,character-method}
\alias{featureInfo<-,ToxicoSet,character,data.frame-method}
\alias{sensitivityRaw,ToxicoSet-method}
\alias{sensitivityRaw<-,ToxicoSet,array-method}
\alias{sensitivitySlot,ToxicoSet-method}
\alias{sensitivitySlot<-,ToxicoSet,list-method}
\alias{sensitivityInfo,ToxicoSet-method}
\alias{sensitivityInfo<-,ToxicoSet,data.frame-method}
\alias{sensitivityProfiles,ToxicoSet-method}
\alias{sensitivityProfiles<-,ToxicoSet,data.frame-method}
\alias{sensitivityProfiles<-,ToxicoSet,matrix-method}
\alias{sensitivityMeasures,ToxicoSet-method}
\alias{drugNames,ToxicoSet-method}
\alias{drugNames<-,ToxicoSet,character-method}
\alias{cellNames,ToxicoSet-method}
\alias{cellNames<-,ToxicoSet,character-method}
\alias{fNames,ToxicoSet,character-method}
\alias{fNames<-,ToxicoSet,character,ANY-method}
\alias{dateCreated,ToxicoSet-method}
\alias{datasetType,ToxicoSet-method}
\alias{datasetType<-,ToxicoSet-method}
\alias{pertNumber,ToxicoSet-method}
\alias{sensNumber,ToxicoSet-method}
\alias{pertNumber<-,ToxicoSet,array-method}
\alias{sensNumber<-,ToxicoSet,matrix-method}
\alias{annotation,ToxicoSet-method}
\alias{annotation<-,ToxicoSet,list-method}
\alias{curation,ToxicoSet-method}
\alias{curation<-,ToxicoSet,list-method}
\title{Class to contain Toxico-genomic Data}
\usage{
\S4method{cellInfo}{ToxicoSet}(object)

\S4method{cellInfo}{ToxicoSet,data.frame}(object) <- value

\S4method{drugInfo}{ToxicoSet}(object)

\S4method{drugInfo}{ToxicoSet,data.frame}(object) <- value

\S4method{phenoInfo}{ToxicoSet,character}(object, mDataType)

\S4method{phenoInfo}{ToxicoSet,character,data.frame}(object, mDataType) <- value

\S4method{molecularProfiles}{ToxicoSet}(object, mDataType, assay)

\S4method{molecularProfiles}{ToxicoSet,character,character,matrix}(object, mDataType, assay) <- value

\S4method{molecularProfiles}{ToxicoSet,character,missing,matrix}(object, mDataType, assay) <- value

\S4method{molecularProfilesSlot}{ToxicoSet}(object)

\S4method{molecularProfilesSlot}{ToxicoSet,ANY}(object) <- value

\S4method{featureInfo}{ToxicoSet,character}(object, mDataType)

\S4method{featureInfo}{ToxicoSet,character,data.frame}(object, mDataType) <- value

\S4method{sensitivityRaw}{ToxicoSet}(object)

\S4method{sensitivityRaw}{ToxicoSet,array}(object) <- value

\S4method{sensitivitySlot}{ToxicoSet}(object)

\S4method{sensitivitySlot}{ToxicoSet,list}(object, ...) <- value

\S4method{sensitivityInfo}{ToxicoSet}(object)

\S4method{sensitivityInfo}{ToxicoSet,data.frame}(object) <- value

\S4method{sensitivityProfiles}{ToxicoSet}(object)

\S4method{sensitivityProfiles}{ToxicoSet,data.frame}(object) <- value

\S4method{sensitivityProfiles}{ToxicoSet,matrix}(object) <- value

\S4method{sensitivityMeasures}{ToxicoSet}(object)

\S4method{drugNames}{ToxicoSet}(object)

\S4method{drugNames}{ToxicoSet,character}(object) <- value

\S4method{cellNames}{ToxicoSet}(object)

\S4method{cellNames}{ToxicoSet,character}(object) <- value

\S4method{fNames}{ToxicoSet,character}(object, mDataType)

\S4method{fNames}{ToxicoSet,character,ANY}(object, mDataType) <- value

\S4method{dateCreated}{ToxicoSet}(object)

\S4method{datasetType}{ToxicoSet}(object)

\S4method{datasetType}{ToxicoSet}(object) <- value

\S4method{pertNumber}{ToxicoSet}(object)

\S4method{sensNumber}{ToxicoSet}(object)

\S4method{pertNumber}{ToxicoSet,array}(object) <- value

\S4method{sensNumber}{ToxicoSet,matrix}(object) <- value

\S4method{annotation}{ToxicoSet}(object)

\S4method{annotation}{ToxicoSet,list}(object) <- value

\S4method{curation}{ToxicoSet}(object)

\S4method{curation}{ToxicoSet,list}(object) <- value
}
\arguments{
\item{object}{A \code{ToxicoSet} object to modify}

\item{value}{A \code{matrix} of replacement values}

\item{mDataType}{\code{character} A string specifying the type of molecular
data to list the phenotype information for.}

\item{assay}{\code{character} Name or index of the assay data to return}

\item{...}{A \code{list} containing fall through arguments; this allows
addition of new parameters to methods for this generic}
}
\value{
An object of the ToxicoSet class

a \code{data.frame} with the cell annotations

Updated \code{ToxicoSet}

a \code{Dframe} with the experiment info

The updated \code{ToxicoSet}

Updated \code{ToxicoSet}

A \code{list} containing the molecularProfiles from a tSet

A copy of the \code{ToxicoSet} with the molecularProfiles slot updated

Updated \code{ToxicoSet}

a \code{data.frame} with the experiment info

Updated \code{ToxicoSet}

a \code{data.frame} with the experiment info

Updated \code{ToxicoSet}

A \code{character} vector of all the available sensitivity measures

A vector of the cell names used in the ToxicoSet

Updated \code{ToxicoSet}

A \code{character} vector of the feature names

Updated \code{ToxicoSet}

The date the ToxicoSet was created

A 3D \code{array} with the number of perturbation experiments per
  radiation type and cell line, and data type

A \code{data.frame} with the number of sensitivity experiments per
  drug and cell line

The updated \code{ToxicoSet}

The updated \code{ToxicoSet}
}
\description{
The ToxicoSet (tSet) class was development to contain and organise large
ToxicGenomic datasets as well as provide useful tools for interacting with
this data. Functions are included for exploring the relationship between
survival fraction and gene expression in cultured human and rat tissues
during exposure to a wide ranges of compounds. Features include plotting dose
and exposure time curves, calculating AUC, fitting linear models and
computing sensitivity signatures.

Get the cell line annotations in a ToxicoSet

Set cell line annotations for a ToxicoSet object

Get the phenotype annotations for cell lines with the specificed molecular
  data type

Set the phenotype annotations for cell lines with the selected molecular
  data type.

Get the molecular profile data associated with the specific molecular data

Set the molecular profile data assciated with the specificied molecular
  data type

Set the feature annotations for a specficied molecular data type

Get the annotations for the sensitivity experiments in the ToxicoSet

Set the annotations for sensitivity experiments in this ToxicSet

Get the data for sensitivty experiments on cell lines in a ToxicoSet

Set the data for sensitivity experiments on cell lines in a ToxicoSet

Get the avilable measurments for sensitivity experiments in a ToxicoSet

Get names of cell lines in a ToxicoSet object

Set the cell line names in a TocicoSet object

Get the feature names in a ToxicoSet object for the specified molecular data
  type

Set the feature names in a ToxicoSet object for the specified molecular
  data type

Get the date a ToxicoSet object was created

Get an array of the number of pertubration experiments per drug and cell
  line in a ToxicoSet object

Get the number of sensitivity experiments per drug and cell line in a
  ToxicoSet

Set the number of perturbation experiments per drug and cell line and
  molecular data type in a ToxicoSet object

Set the number of sensitivity experiments per drug and cell line in a
  ToxicoSet object
}
\section{Methods (by generic)}{
\itemize{
\item \code{cellInfo}: Returns the annotations for all the cell lines tested
on in the ToxicoSet

\item \code{cellInfo<-}: Returns the annotations for all the cell lines tested
on in the ToxicoSet

\item \code{drugInfo}: Returns the annotations for all the drugs tested in
the ToxicoSet

\item \code{drugInfo<-}: Update the drug annotations

\item \code{phenoInfo}: Return the experiment info from the given type of
molecular data in ToxicoSet

\item \code{phenoInfo<-}: Update the the given type of molecular data experiment
info in the ToxicoSet

\item \code{molecularProfiles}: Return the given type of molecular data from the
ToxicoSet

\item \code{molecularProfiles<-}: Update the given type of molecular data from the
ToxicoSet

\item \code{molecularProfiles<-}: Update the given type of molecular data from the
ToxicoSet

\item \code{molecularProfilesSlot}: Get contents of molecularProfiles slot

\item \code{molecularProfilesSlot<-}: Update the molecular profiles slot of a ToxicoSet and
returns the updated copy

\item \code{featureInfo}: Return the feature info for the given molecular data

\item \code{featureInfo<-}: Replace the gene info for the molecular data

\item \code{sensitivityRaw}: Retrive the raw dose and viability data from an tSet

\item \code{sensitivityRaw<-}: Set the raw dose and viability data for a tSet and return
and updated copty

\item \code{sensitivitySlot}: Retrieves the contents of the sensitivity slot

\item \code{sensitivitySlot<-}: Set the raw dose and viability data for an tSet and return
and updated copy

\item \code{sensitivityInfo}: Return the drug dose sensitivity experiment info

\item \code{sensitivityInfo<-}: Update the sensitivity experiment info

\item \code{sensitivityProfiles}: Return the phenotypic data for the drug dose
sensitivity

\item \code{sensitivityProfiles<-}: Update the phenotypic data for the drug dose
sensitivity

\item \code{sensitivityProfiles<-}: Update the phenotypic data for the drug dose
sensitivity

\item \code{sensitivityMeasures}: Returns the available sensitivity profile
summaries, for example, whether there are IC50 values available

\item \code{drugNames}: Return the names of the drugs used in the ToxicoSet

\item \code{drugNames<-}: Update the drug names used in the dataset

\item \code{cellNames}: Return the cell names used in the dataset

\item \code{cellNames<-}: Update the cell names used in the dataset

\item \code{fNames}: Return the feature names used in the dataset

\item \code{fNames<-}: Update the feature names used in the dataset

\item \code{dateCreated}: Return the date the ToxicoSet was created

\item \code{datasetType}: Update the dataset type of an tSet and return a copy of
the updated object

\item \code{datasetType<-}: Update the dataset type of an tSet and return a copy of
the updated object

\item \code{pertNumber}: Return the summary of available perturbation
experiments

\item \code{sensNumber}: Return the summary of available sensitivity
experiments

\item \code{pertNumber<-}: Update the summary of available perturbation
experiments

\item \code{sensNumber<-}: Update the summary of available sensitivity
experiments

\item \code{annotation}: Retrieve the annotations slot form an tSet

\item \code{annotation<-}: Update the annotation slot of a tSet

\item \code{curation}: Retrieve the curation slot form an tSet

\item \code{curation<-}: Update the annotation slot of a tSet
}}

\section{Slots}{

\describe{
\item{\code{annotation}}{A \code{list} of annotation data about the ToxicoSet,
including the \code{$name} and the session information for how the object
was creating, detailing the exact versions of R and all the packages used}

\item{\code{molecularProfiles}}{A \code{list} containing \code{SummarizedExperiment}
type object for holding data for RNA, DNA, SNP and CNV
measurements, with associated \code{fData} and \code{pData}
containing the row and column metadata}

\item{\code{cell}}{A \code{data.frame} containing the annotations for all the cell
lines profiled in the data set, across all data types}

\item{\code{drug}}{A \code{data.frame} containg the annotations for all the drugs
profiled in the data set, across all data types}

\item{\code{sensitivity}}{A \code{list} containing all the data for the sensitivity
experiments, including \code{$info}, a \code{data.frame} containing the
experimental info,\code{$raw} a 3D \code{array} containing raw data,
\code{$profiles}, a \code{data.frame} containing sensitivity profiles
statistics, and \code{$n}, a \code{data.frame} detailing the number of
experiments for each cell-drug pair}

\item{\code{perturbation}}{A \code{list} containting \code{$n}, a \code{data.frame}
summarizing the available perturbation data,}

\item{\code{curation}}{A \code{list} containing mappings for \code{$drug},
\code{cell}, \code{tissue} names  used in the data set to universal
identifiers used between different ToxicoSet objects}

\item{\code{datasetType}}{A \code{character} string of 'sensitivity',
'perturbation', or both detailing what type of data can be found in the
ToxicoSet, for proper processing of the data}
}}

\examples{
data(TGGATESsmall)
cellInfo <- cellInfo(TGGATESsmall)

data(TGGATESsmall)
cellInfo(TGGATESsmall) <- cellInfo(TGGATESsmall)

data(TGGATESsmall)
phenoInfo <- phenoInfo(TGGATESsmall, mDataType="rna")

data(TGGATESsmall)
phenoInfo(TGGATESsmall, mDataType="rna") <-
  phenoInfo(TGGATESsmall, mDataType="rna")

data(TGGATESsmall)
TGGATES_mProf <- molecularProfiles(TGGATESsmall, "rna")[seq_len(10),]

molecularProfiles(TGGATESsmall, "rna") <-
  molecularProfiles(TGGATESsmall, "rna")

data(TGGATESsmall)
molecularProfilesSlot(TGGATESsmall)

data(TGGATESsmall)
molecularProfilesSlot(TGGATESsmall) <- molecularProfilesSlot(TGGATESsmall)

data(TGGATESsmall)
featureInfo <- featureInfo(TGGATESsmall, "rna")[seq_len(10),]

data(TGGATESsmall)
featureInfo(TGGATESsmall, "rna") <- featureInfo(TGGATESsmall, "rna")

data(TGGATESsmall)
sensitivitySlot(TGGATESsmall)

data(TGGATESsmall)
sensitivitySlot(TGGATESsmall) <- sensitivitySlot(TGGATESsmall)

data(TGGATESsmall)
sensInf<- sensitivityInfo(TGGATESsmall)[seq_len(10),]

data(TGGATESsmall)
sensitivityInfo(TGGATESsmall) <- sensitivityInfo(TGGATESsmall)

data(TGGATESsmall)
sensProf <- sensitivityProfiles(TGGATESsmall)

sensitivityProfiles(TGGATESsmall) <- sensitivityProfiles(TGGATESsmall)

sensitivityMeasures(TGGATESsmall)

cellNames(TGGATESsmall)

data(TGGATESsmall)
cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)

fNames(TGGATESsmall, "rna")[seq_len(10)]

data(TGGATESsmall)
cellNames(TGGATESsmall) <- cellNames(TGGATESsmall)

dateCreated(TGGATESsmall)

data(TGGATESsmall)
datasetType(TGGATESsmall)

pertNumber(TGGATESsmall)

sensNumber(TGGATESsmall)

pertNumber(TGGATESsmall) <- pertNumber(TGGATESsmall)

sensNumber(TGGATESsmall) <- sensNumber(TGGATESsmall)

}
