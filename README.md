
# CensoredMediator

<!-- badges: start -->
<!-- badges: end -->

The goal of CensoredMediator is analyze the mediation model with a
censored mediator, described in:

“Wang J, Ning J, Shete S. Mediation analysis in a case-control study
when the mediator is a censored variable. Stat Med 38(7):1213-1229,
3/2019. PMCID: PMC6467083”

“Wang J, Ning J, Shete S. Mediation model with a categorical exposure
and a censored mediator with application to a genetic study. PLoS One
16(10):e0257628, 2021. e-Pub 10/2021. PMCID: PMC8509986”

## Installation

You can install the development version of CensoredMediator like so:

``` r
library(devtools)
devtools::install_github("lilyxj91/CensoredMediator",force = TRUE)
```

## Example

``` r
library(CensoredMediator)
#> Warning: replacing previous import 'MASS::select' by 'dplyr::select' when
#> loading 'CensoredMediator'
## basic example code
```

### Read the data

``` r
### Read the data
data<-TestData
head(data)
#>   Y X        M delta          Z1 Z2
#> 1 0 1 5.865867     1 -0.09050880  0
#> 2 0 0 6.873086     1  0.32238550  1
#> 3 0 1 2.944863     0 -0.20819192  0
#> 4 0 1 6.288533     1 -0.01466635  1
#> 5 1 0 6.139924     1  0.22106215  1
#> 6 0 0 5.802920     1 -0.14373902  1

# where 
# Y is the outcome;
# X is the initial variable;
# M is the censored mediator;
# delta is the censored status for M (1-observed; 0-censored);
# Z1 and Z2 are two covariates.  
# 
# Please arrange your data into this format, with the first four columns named as Y, X, M and delta; and then include other covariates from column five (Z1, Z2, Z3...).  
```

### Define the parameters

``` r
## If the outcome is a continuous variable: 
## use the following parameters
# outcome.type<-"continuous"
# casecontrol<-FALSE
# prev<-0
# boot<-200  ## number of bootstrap samples
```

``` r
## If the outcome is a binary variable but the study is not a case-control study: 
## use the following parameters
outcome.type<-"binary"
casecontrol<-FALSE
prev<-0
boot<-200  ## number of bootstrap samples
```

``` r
## If the outcome is a binary variable and the study is a case-control study:
## use the following parameters
# outcome.type<-"binary"
# casecontrol<-TRUE
# prev<-0.12 ## prevalence of the disease (case) in general population; non-zero for a case-control study
# boot<-200  ## number of bootstrap samples
```

### Calculate the coefficients for each paths

``` r
(results.tmp<-CensoredMediator(data,outcome.type,casecontrol,prev,printout=TRUE))
#> Coefficients of M~X:
#> 
#> (Intercept)           X          Z1          Z2 
#>       5.993       0.393       0.318       0.263 
#> 
#> 
#> Coefficients of Y~X+M:
#> 
#> (Intercept)           X           M          Z1          Z2 
#>      -5.182       0.808       0.373       0.610       0.545
#> [[1]]
#>      (Intercept)         X        Z1        Z2
#> [1,]    5.992783 0.3932519 0.3178537 0.2629342
#> 
#> [[2]]
#> (Intercept)           X           M          Z1          Z2 
#>  -5.1821041   0.8084458   0.3734420   0.6099168   0.5449263 
#> 
#> [[3]]
#>           IE         DE         TE        PM
#> 1 0.01925094 0.07716508 0.09641602 0.1996654

(IE.results<-results.tmp[[3]])
#>           IE         DE         TE        PM
#> 1 0.01925094 0.07716508 0.09641602 0.1996654

### Calculate the indirect effect (IE) and percentage mediated (PM)
### and their corresponding 95% confidence intervals using bootstrapping
### May take hours or days to run depending on the number of bootstrap samples
(results.boot<-mybootstrap(data,IE.results,boot,outcome.type,casecontrol,prev,printout=FALSE))
#> Indirect effect (IE) and 95% CI:
#> 
#> [1] "0.019 (0.008, 0.031)"



#> Percentage mediated (PM) and 95% CI: 
#> 
#> [1] "0.2 (0.091, 0.358)"
#> $CI_IE
#> [1] 0.008098255 0.031216266
#> 
#> $CI_PM
#> [1] 0.09128407 0.35818356
