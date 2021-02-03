# Analysis code for "CD4 receptor diversity represents an ancient protection mechanism against primate lentiviruses" 

To regenerate the analysis, run `source('virusAlleleAnalysis.R')` in R.

The code depends on the R package `rstan` available from CRAN (`install.packages('rstan')`) and was run in R version 3.6.3 with `rstan` version 2.21.1.

This code should generate analysis result files in the `out` directory. Example output is available in the `exampleOut` directory where:

* BonoboDiffs.csv: contains Bayesian posterior estimates for comparisons of the infectivity of the various Env in bonobo alleles. Columns are:
  1. Virus Virus from which the Env was derived
  2. AlleleA First allele in the comparison
  3. AlleleB Second allele in the comparison
  4. p(fold change < 1): The posterior probability that the fold difference in odds between allele A and allele B is less than 1
  5. Estimated fold change: The posterior mean estimate for the fold change in odds between allele A and allele B
  6. Lower 95% CrI: The lower limit of the estimated 95% credible interval for the fold change in odds between allele A and allele B
  7. Upper 95% CrI: The upper limit of the estimated 95% credible interval for the fold change in odds between allele A and allele B
* OtherDiffs.csv: As above but for other primate alleles
* GorillaDiffs.csv: As above but for gorilla alleles
* N15T.csv: As above but comparing the estimated effects of changing N15T in the AHSM gorilla allele


