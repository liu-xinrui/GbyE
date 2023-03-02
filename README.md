# **GbyE**
Effect model of gene-environment interaction, involved in predictive analysis of GWAS and GS<br>
# **Introduction**
Gbye main package, coordinate other function packages for GWAS and GS operations, sometimes, We need to prepare documents, including <br>
<br>
# **Improve GbyE's test function package**<br>
Here is a small function package related to gbye operation, which helps to run gbye programs more conveniently<br>
This currently includes:<br>
### `GbyE.R`
Gbye main package, coordinate other function packages for GWAS and GS operations, sometimes, We need to prepare documents, including <br>
   * `GD` (Genotype file,type [0,1,2])
   * `GM` (SNP marker map file with SNP names, Chrosome and Position)
   * `Y` (Phenotype file,If the model phenotype is not used and the real phenotype prediction is directly used, please enter this file)
   * `ha2` (Heritability of the simulate phenotype)
   * `cov_g` (Genetic correlation of the simulate phenotyp)
   * `NQTN` (Number of phenotype related SNPs simulated, defult=20)
   * `nrep` (Number of iterations)
   * `nfold` (CV Multiple cross validation)
   * `nIter` & `burnIn` (set Markov chain length of the R packages "**BGLR**")
   * `gwas` (defult=T, Whether to perform GWAS operation)
   * `gs` (defult=T, Whether to perform GS operation)
   * `file.output` (Whether the operation result file is output)
   * `plot` (Whether curve drawing is carried out on the results can show the prediction accuracy of the results)
### `G&E_Simulation.R`
Function package for phenotypic simulation according to genotype file<br>
### `GbyE.file.Calculate.R`
Gbye file is generated from the original genotype file and SNP map file to show the effect of gene environment interaction<br>
### `Power.FDR.Calculate.R`
Calculate the FDR and Power values of the prediction results<br>
### `Comparison.GWAS.Result.PValue.R`
In GWAS operation, compare the pvalue value of the same SNP and screen the most significant p value for subsequent GS prediction<br>
