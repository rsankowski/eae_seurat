# adjusted from Josip Herman
library(readr)
library(data.table)
library(tidyverse)
library(tools)
library(assertthat)

#this code was adjusted from Josip Herman, MPI Freiburg
# Specify directories
file_paths <- list(file.path("data","counts"))
names(file_paths) <- file_paths


# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".coutt.csv")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  file.path(names(files_by_ext[x]), files_by_ext[[x]]) } ))

names(all_file_paths) <- gsub("data/counts/|.coutt.csv", "", all_file_paths)

# Calculate md5sums to check for duplicated
md5sums      <- lapply(all_file_paths, function(x) {md5sum(x)} )

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1))))) == 0)

####
#### LOADING
####
# Loading data using lapply
data_list   <- lapply(all_file_paths, function(x) {fread(x, header= T)} )

# Add dataset name prefix to all columns, Merge with remaining gene names
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))
  
}

# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- reduce(data_list, full_join, by = "GENEID")
data_list_cbind[is.na(data_list_cbind)]   <- 0

# Remove ERCC and mitochondrial genes
prdata <- as.data.frame(data_list_cbind)
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL

#remove Kcnq1ot1 and correlated genes
cs <- colSums(prdata)
prdata <- prdata[,cs > 499]

#exclude not detected genes
prdata <- prdata[rowSums(prdata) >0,]

cs <- cs[cs > 499]
f <- t(prdata["Kcnq1ot1__chr11",])/cs < .02

prdata <- prdata[!duplicated(gsub("_.*", "", rownames(prdata))),]
rownames(prdata) <- gsub("_.*", "", rownames(prdata))

save(prdata, file = "data/prdata.RData")
