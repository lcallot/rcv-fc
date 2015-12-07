# Computation material for: 'Modeling and forecasting large realized covariance matrices and portfolio choice.'
## Laurent Callot, Anders B. Kock, Marcelo C. Medeiros.


[Link to the paper](http://lcallot.github.io/papers/rcv_forecasts.pdf/)

[Link to the supplementary material](http://lcallot.github.io/papers/rcv_supplement.pdf/)

---
Date 06/12/2015



This repository contains the material used to compute the results in: 'Modeling and forecasting large realized covariance matrices and portfolio choice.' Due to the complexity of the computations involved and the size of the intermediate output, the material in the present repo does not contain the saved forecast results necessary to generate the output. All the code used in the paper as well as the raw data is included.


### Content

Our computations can be summarized in 3 main steps. 

1. Computing the VAR forecasts. This step involves heavy computations, some of which was carried on the [Lisa cluster](https://userinfo.surfsara.nl/systems/lisa) and some on the server of the VU's Econometrics and OR department. The **R** scripts for these computations are found in the _VAR\_forecasting_ folder. 
2. Computing summary statistics of the forecasts and generating tables and figures. The **knitr** (.Rnw) files used for this step can be found in the _VAR\_processing_ folder. Generating the output requires that the results from step 1 be included in the _fc\_data_ folder. The output files of the forecasting procedure are not included in this repo due to their large size. 
3. The files used to compute the portfolio statistics. These files are cound in the _portfolio_ folder, the main file is _portopt.m_. The input for the portfolio selection procedure is _.csv_ files containing the stacked forecasts of the covariance matrix and files containing the daily returns found in the _data_ folder. The output is saved in the _.mat_ file.  Note that these are **Matlab** files written by Marcelo C. Medeiros, they require the _optimization toolbox_.

Other material includes:

 - __data__: The realized covariance matrices, daily returns, dates of the observations, and names and industry categories for the stocks of the Dow Jones. 
 - __subs__: Sets of __R__ functions used in the computations.
 - __DCC\_and\_EWMA__: The code used to compute the DCC and EWMA models. 





### Required packages 



```r
library('devtools')
install_github('lcallot/lassovar')
library('lassovar')

library(plyr)
library(parallel)
library(magrittr)
library('Matrix')
library('SparseM')
library('glmnet')
library('xtable')
library('expm')
library('doMC')

# for plotting
library('reshape2')
library('ggplot2')

# for the DCC
library(rmgarch)

```

All packages except __lassovar__ are available from CRAN. 



