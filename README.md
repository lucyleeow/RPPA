# RPPA

R package for performing basic data cleaning, exploratory visualisations and normalisation of reverse phase protein array (RPPA) data from the [Zeptosens system](https://www.ncbi.nlm.nih.gov/pubmed/12164697). RRPA is an antibody based assay used to determine protein expression levels.

The output from the the Zeptosens system are relative fluorescence intensity (RFI) values, which are indicative of protein concentration in that sample, for each antibody and protein. The format of the data is a matrix where each row represents a sample and each column represents a antibody.

To intall this package, first make sure you have installed and loaded [devtools](https://cran.r-project.org/web/packages/devtools/readme/README.html) (using `install.packages("devtools")` then `library(devtools)`). 

You can then install using:
`install_github("lucyleeow/RPPA")`
and load with:
`library(RPPA)`

