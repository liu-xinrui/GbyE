# **GbyE**
Effect model of gene-environment interaction, involved in predictive analysis of GWAS and GS<br>
![GbyE](https://raw.githubusercontent.com/liu-xinrui/GbyE/main/base/GbyE.png)
## **Introduction**
Gbye main package, coordinate other function packages for GWAS and GS operations, sometimes, We need to prepare documents, including <br>
   * `GD` (Genotype file,type [0,1,2])
   * `GM` (SNP marker map file with SNP names, Chrosome and Position)
   * `Y` (Phenotype file,The phenotype file must contain two environments or two traits)
   * `GbyE.GD` & `GbyE.GM` (If you have previously performed GbyE operations and saved a GbyE file, you can add the file directly to GbyE to reduce the time cost incurred by calculating again)
   * `gwas` (defult=T, Whether to perform GWAS operation)
   * `PCA.total` (Number of principal components as covariates at GWAS)
   * `CV` (GbyE is based on GAPIT operations, and like GAPIT, it runs with covariates for operations)
   * `gwas.model`(Select the GWAS analysis model you need, theoretically you can use all GAPIT models, such as: '___MLM___', '___CMLM___', '___BLINK___', '___MLMM___')
   * `cutOff`(Allows setting a threshold of significance when drawing)
   * `gs` (defult=T, Whether to perform GS operation)
   * `method`(Set the operation method of GS, currently only '___BGLR___', '___rrBLUP___' and '___GAPIT___' are supported)
   * `gs.model`(GS model: rrblup-default-"___ML___", gapit "___gblup___" "___mablup___" "___cblup___" "___sblup___", bglr "___RRB___" "___BA___" "___BC___" "___BB___" "___BL___")
   * `nIter` & `burnIn` (set Markov chain length of the R packages "**BGLR**")
   * `file.output` (Whether the operation result file is output)
   * `plot` (Whether curve drawing is carried out on the results can show the prediction accuracy of the results)
   * `file.type`(Set the output image format, the default is PDF)
   * `dpi`(Set the resolution of the output image, default=600)
## **Output**
The output results contain three panels, `GbyE`, `GWAS` and `GS`. `GbyE` mainly includes the output GbyE genotype file and phenotype file. `GWAS` mainly includes additive effect `adde`, reciprocal effect `GXE`, overall `summary` result of GWAS, `GMM` mapping parameter file. `GS` contains mainly, prediction results of all missing values `pred`. The prediction parameters of `rrBLUP`. The prediction parameters of GAPIT software, `gapit`. The prediction parameters of BGLR software, `bglr`.
![result tree](https://raw.githubusercontent.com/liu-xinrui/GbyE/main/base/GbyE.result.tree.jpg)
# **Predictable results**
Gbye main package, coordinate other function packages for GWAS and GS operations, sometimes, We need to prepare documents, including <br>
## **Genome-wide association studies**
***Manhattan Plots*** <br>
![Manhattan Plot by demo data](https://raw.githubusercontent.com/liu-xinrui/GbyE/main/base/Manhattan.jpg)
***Quantile-Quantile Plots*** <br>
![QQ Plot by demo data](https://raw.githubusercontent.com/liu-xinrui/GbyE/main/base/QQplot.jpg)
<br>
# **How to use**
```
#Loading support packages
source("https://raw.githubusercontent.com/liu-xinrui/GbyE/main/GbyE.R")
source("https://raw.githubusercontent.com/liu-xinrui/data/main/gapit_functions.txt")
source("https://raw.githubusercontent.com/liu-xinrui/data/main/GAPIT.library.R")

#Please import your own data below
GD=read.table("mdp_numeric.txt",head=T)
GM=read.table("mdp_SNP_information.txt",head=T)
Y=read.table("mdp_traits.txt",head=T)

#Remove phenotypic data without genotypes to maintain data consistency
Y=Y[match(GD[,1],Y[,1]),]

#Run GbyE
myGbyE=GbyE(GD=GD,
            GM=GM,Y=Y,
            PCA.total=3,
            gwas=F,gs=T,
            plot=T,
            gwas.model="MLM",
            method="gapit")
```

# **Improve GbyE's test function package**<br>
Here is a small function package related to gbye operation, which helps to run gbye programs more conveniently<br>
This currently includes:<br>
### `GbyE.R`
The following is used to support GbyE runs and to run tests and refine GbyE builds. unlike the main GbyE function package, GbyE.test.R adds simulation functions for testing GbyE only. the following are the parameters of this function, including <br>
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
