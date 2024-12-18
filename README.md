# QTL2-implementation
- A shallow implementation of Karl Broman's QTL2 R package for Quantative Trait Loci (QTL) analysis.
- QTL2 package available at https://github.com/kbroman/qtl2
- This script is intended for use with Mus Musculus; however, you could always just change the BioMart mart.

## The results
- A PDF that shows QTL peaks by tested phenotype. Test alphas are set to .20 (yellow-green), .10 (orange) , and .05 (red) by default. 
- A CSV that is comprised of queried Ensembl data within significant QTL peaks.
- A CSV that shows the numbers that form the basis of the QTL peaks graphed on PDF.

## How to use?
- Simply replace the default value for the .yaml control value with the actual file name of your control file
- Change any file names, passed as arguments to custom functions, to personalize the results.

## Use this reference control file to guide the creation of your control file.
 - https://github.com/kbroman/qtl2/blob/gh-pages/assets/sampledata/iron/iron.yaml
   
## Screenshots
![sanitized](https://github.com/user-attachments/assets/5020cde9-dd16-47cc-a2d7-e0635f4a8c83)
