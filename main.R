library(tidyverse)

source("R/PhosphositeAnalysis_c2a.R")
source("R/KinaseAnalysis.R")

if (!dir.exists("results")){
  dir.create("results")
}

# 1. get Phosphosite files

read_psites_from_folder(folder = data1, foldernum = '1')
read_psites_from_folder(folder = data2, foldernum = '2')

# read_psites_from_folder(folder = "../compare2analyses_data/example1", foldernum = '1')
# read_psites_from_folder(folder = "../compare2analyses_data/example2", foldernum = '2')

# Data checks
if (exists('lfcs1') & exists('lfcs2')){
  # Check if each BN LFC file has a matching p file in a folder
  print("Checking the names of BN files...")
  are_lfcs_ps_paired(lfcs1, ps1, foldernum = '1')
  are_lfcs_ps_paired(lfcs2, ps2, foldernum = '2')
  
  # Check if each BN LFC and p file has a matching LFC and p file in the other folder
  is_data_in_2folders_paired_bn(lfcs1, lfcs2, ps1, ps2)
  print("Names of BN files are matching in both folders.")
  
}

if (exists('tercen1') & exists('tercen2')){
  is_data_in_2folders_paired_tercen(tercen1, tercen2)
  print("Names of Tercen files are matching in both folders.")
}

# 2. get peptide results

if (exists('lfcs1')){
  get_peptide_results_bn(lfcs1, ps1, lfcs2, ps2, pcutoff = pcutoff,
                         data1_name = data1_name, data2_name = data2_name)
}

if (exists('tercen1')){
  get_peptide_results_tercen(tercen1, tercen2, pcutoff = pcutoff,
                             data1_name = data1_name, data2_name = data2_name)
}

# 3. get UKA files

uka_files1 <- list.files(path = data1, pattern = "UKA_", full.names = TRUE)
uka_files2 <- list.files(path = data2, pattern = "UKA_", full.names = TRUE)

# uka_files1 <- list.files(path = "../compare2analyses_data/example1", pattern = "UKA_", full.names = TRUE)
# uka_files2 <- list.files(path = "../compare2analyses_data/example1", pattern = "UKA_", full.names = TRUE)


stopifnot("Number of UKA files in data1 and data2 folders don't match!" =
            length(uka_files1) == length(uka_files2))

# 4. Get UKA results

get_uka_results(uka_files1 = uka_files1, uka_files2 = uka_files2, output_type = output_type, data1_name = data1_name, 
                data2_name = data2_name, fscorecutoff = fscorecutoff, uka_version = uka_version)







