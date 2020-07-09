# R dag 3

```R
library(MESS)
library(glmnet)
mres <- mfastLmCpp(pheno$V3, as.matrix(geno2))
summary(phenotypes)
mres <- mfastLmCpp(phenotypes$BMI,   (genotypes))
summary(mres)
pval <- 2*pt(-abs(mres[,3]), df=nrow(phenotypes)-2)
summary(pval)
head(sort(pval))
head(names(genotypes)[order(pval)]) #Hvis genotypes no names, remove "names".
summary(phenotypes)
mres <- mfastLmCpp(phenotypes$BMI, as.matrix(genotypes))
pval <- 2*pt(-abs(mres[,3]), df=nrow(phenotypes)-2)
summary(mres)
length(mres)
length(mres$coefficients)
length(phenotypes$BMI)
statistics <- mres[,3]^2 ; median(statistics)
median(statistics)
rescaled <- statistics/(median(statistics)/0.456)
head(rescaled)
results <- 1-pchisq(rescaled, df=1)
head(sort(results))
plot(-log(results))
res <- glmnet(genotypes, phenotypes$BMI)
plot(res) # find lambda sd1 og smid i glmnet
# manhattan kan ikke plottes fra lasso
head(mres)
  coefficients     stderr     tstat
1  -0.16280396 0.07554746 -2.154989
2  -0.08726746 0.06026324 -1.448104
3  -0.08526341 0.07970901 -1.069683
4  -0.06899764 0.05790543 -1.191557
5  -0.08779398 0.07742912 -1.133863
6  -0.01763163 0.05534878 -0.318555
> head(sort(pval))
[1] 3.502032e-07 1.301984e-06 4.814747e-06 1.202277e-05 1.268807e-05
[6] 1.901674e-05
> plot(pvals
+      )
Error in plot(pvals) : object 'pvals' not found
> plot(pval
+      )
> plot(-log(pval))
> log(1)
[1] 0
> log10(1)
[1] 0
> chromosome <- rep(1:2, times=c(11283,32019-11283))
> summary(chromosome)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0     1.0     1.5     1.5     2.0     2.0 
> plot(-log(pval),col=chromosome)
> plot(-log(results),col=chromosome)
> plot(-log(results),col=chromosome)

head(order(pval)) #= SNPs of interest

> head(sort(statistic1))
[1] -134.24022  -84.83124  -79.05183  -78.17256  -56.61404  -56.17081
> statistic1 <- mres[,3]^2
> head(sort(statistic1))
[1] 7.281995e-09 2.171639e-08 2.872567e-08 3.185324e-08 3.193815e-08 4.256086e-08
> median(statistic1)
[1] 0.4544575
> plot(-log(pval),-log(seq(0,1,length=32021)[-c(1,32021)]))
> plot(sort(-log(pval)),-log(sort(seq(0,1,length=32021)[-c(1,32021)])))
> plot(sort(-log(pval)),sort(-log(seq(0,1,length=32021)[-c(1,32021)])))
> abline(0,1,col='red')
```

Informationen man får fortæller kun p-værdi som er afhængig af n, siger intet om klinisk relevans. Effect size er vigtigt. Kan ikke ses på Manhattan plot.



Change close to median, needs correction, do you believe that 40% of genes affect the given outcome.



glmnet bruger mere effect size. 

pcc etc. = marginal analysis.

Flere genotypes der ligger tæt fjernes af lasso

