# General
Compares peptide and kinase results of two projects with the same design (repeated projects).

# File structure
Save files of the two projects in separate folders.
In each folder, the files should have the exact same naming. 
In each file, the namespaces should be exactly the same.

Name files as for automated reporting. 
Examples:
MTvC_PTK_01_TvsC_LogFC.txt
TT_STK_01_TvsC_p.txt
UKA_STK_01_TvsC.txt


# RUN:
Navigate to the folder of the script, where the two data folders to be compared are saved.

Parameters of the run command:
Location of R
name of the R script (compare2analysis.R)
location of folder 1
location of folder 2
name of data 1 (plotted in the x axis), e.g. "July data"
name of data 2 (plotted in the y axis), e.g. "August data"
pcutoff: p-value threshold, default = 0.05
fscorecutoff: Median Final Score threshold, default = 1.3


## Example:
"C:\Program Files\R\R-4.2.3\bin\Rscript.exe" compare2analysis.R --folder1 "./data_ex1/" --folder2 "./data_ex2/" --data1_name "Project1" --data2_name "Project2"



# Features to be developed:

* N_common kinases: like n common peps output table
* N_common peptides csv: colnames from data1, data2 in opt
* UKA Scores: get median instead of mean when available
