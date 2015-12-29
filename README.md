# Computation material for: 'Modeling and forecasting large realized covariance matrices and portfolio choice.'
## Laurent Callot, Anders B. Kock, Marcelo C. Medeiros.


[Link to the paper](http://lcallot.github.io/papers/rcv_forecasts.pdf)

[Link to the supplementary material](http://lcallot.github.io/papers/rcv_supplement.pdf)

---
Date 29/12/2015



This repository contains the material used to compute the results in: 'Modeling and forecasting large realized covariance matrices and portfolio choice.' Due to the complexity of the computations involved and the size of the intermediate output, the material in the present repo does not contain the saved forecast results necessary to generate the output. All the code used in the paper as well as the raw data is included.


### Data 

The data is contained in the _data_ folder. It contains 13 files

  - Files starting with __CRK__ are __R__ data files containing realized covariance matrices for different subset of the data at different levels of aggregation. These are __Rdata__ files and can be loaded with the command `load('CRK_file_name)` in __R__. The file names should be interepreted as follows:
    1. Files names ending with _W_ or _M_ refer to weekly or monthly aggregated data, all others (with _har_ in the name) are daily data. 
    2. Files containing _dj_ in the name indicate that the data is composed of the 30 stocks of the Dow-Jones. Files containing _fac_ are the Dow-Jones augmented with the S\&P 500 used as a common factor. 
    3. Files containing _cens_ in the name refer to censored data (see paper), _none_ is used for uncensored data. 
    4. The transformations applied to the data are refered to as _lcov_ and _lmat_ or _none_, see the paper for details. 
    
 - __mk\_aggdata.Rnw__ and __mk\_aggdata.pdf__ contains the code used to aggregate the daily data to weekly and monthly data.
 - __dates__ and __dj\_crk\_dates.txt__ are plain text files containing the calendar date of the daily observations. 
 - __sp\_indus.csv__ is a plain text file containing information on the industry category of every stock in the S\&P 500. __dj-cat__ and __dj\_crk\_names.txt__ are subsets of that file for the Dow Jones stock. __get\_dj\_indus.R__ is the script used to extract the Dow Jones subset. 
 - __dj-ind__ is a plain text file containing the ticker and S\&P index of the 30 stocks of the Dow Jones. 
 

### Computations

Our computations can be summarized in 3 main steps. 

1. Computing the VAR forecasts. This step involves heavy computations, some of which was carried on the [Lisa cluster](https://userinfo.surfsara.nl/systems/lisa) and some on the server of the VU's Econometrics and OR department. The **R** scripts for these computations are found in the _VAR\_forecasting_ folder. 
2. Computing summary statistics of the forecasts and generating tables and figures. The **knitr** (.Rnw) files used for this step can be found in the _VAR\_processing_ folder. Generating the output requires that the results from step 1 be included in the _fc\_data_ folder. The output files of the forecasting procedure are not included in this repo due to their large size. 
3. The files used to compute the portfolio statistics. These files are cound in the _portfolio_ folder, the main file is _portopt.m_. The input for the portfolio selection procedure is _.csv_ files containing the stacked forecasts of the covariance matrix and files containing the daily returns found in the _data_ folder. The output is saved in the _.mat_ file.  Note that these are **Matlab** files written by Marcelo C. Medeiros, they require the _optimization toolbox_.

Other material includes:

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



