
library(tidyverse)
library(optparse)

source("R/PhosphositeAnalysis_c2a.R")
source("R/KinaseAnalysis.R")


defaultW <- getOption("warn")
options(warn = -1)


option_list = list(
  make_option(c("--folder1"), action="store", default=NA, type='character',
              help="path to data1 folder"),
  make_option(c("--folder2"), action="store", default=NA, type='character',
              help="path to data2 folder"),
  make_option(c("--data1_name"), action="store", default=NA, type='character',
              help="name of earlier project"),
  make_option(c("--data2_name"), action="store", default=NA, type='character',
              help="name of later project"),
  make_option(c("--pcutoff"), action="store", default=0.05, type='numeric',
              help="significance cutoff for peptide data"),
  make_option(c("--fscorecutoff"), action="store", default=1.3, type='numeric',
              help="Median Final score cutoff for kinase data")

)

opt = parse_args(OptionParser(option_list=option_list))

print(opt)


# 1. get Phosphosite files

read_psites_from_folder(folder = opt$folder1, foldernum = '1')
read_psites_from_folder(folder = opt$folder1, foldernum = '2')

read_psites_from_folder(folder = "../compare2analyses_data/example1", foldernum = '1')
read_psites_from_folder(folder = "../compare2analyses_data/example2", foldernum = '2')


# Check if each BN LFC file has a matching p file in a folder
print("Checking the filenames...")
are_lfcs_ps_paired(lfcs1, ps1, foldernum = '1')
are_lfcs_ps_paired(lfcs2, ps2, foldernum = '2')

# Check if each BN LFC and p file has a matching LFC and p file in the other folder
is_data_in_2folders_paired_bn(lfcs1, lfcs2, ps1, ps2)
is_data_in_2folders_paired_tercen(tercen1, tercen2)
print("Filenames are correct.")


# 2. get peptide results
source("R/PhosphositeAnalysis_c2a.R")
get_peptide_results_bn(lfcs1, ps1, lfcs2, ps2, pcutoff = 0.05,
                       data1_name = "data1", data2_name = "data2")

get_peptide_results_tercen(tercen1, tercen2, pcutoff = 0.05,
                       data1_name = "data1", data2_name = "data2")


# 3. get UKA files

uka_files1 <- list.files(path = opt$folder1, pattern = "UKA_", full.names = TRUE)
uka_files2 <- list.files(path = opt$folder2, pattern = "UKA_", full.names = TRUE)

uka_files1 <- list.files(path = "../compare2analyses_data/example1", pattern = "UKA_", full.names = TRUE)
uka_files2 <- list.files(path = "../compare2analyses_data/example1", pattern = "UKA_", full.names = TRUE)


stopifnot("Number of UKA files in data1 and data2 folders don't match!" =
            length(uka_files1) == length(uka_files2))

# 4. Get UKA results
source("R/KinaseAnalysis.R")
get_uka_results(uka_files1, uka_files2, data1_name ="data1", data2_name = 'data2', fscorecutoff)







