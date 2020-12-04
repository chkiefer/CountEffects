# CountEffects
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/chkiefer/CountEffects.svg?branch=master)](https://travis-ci.com/chkiefer/CountEffects)
[![HitCount](http://hits.dwyl.com/chkiefer/CountEffects.svg)](http://hits.dwyl.com/chkiefer/CountEffects)

R package to computed average and conditional effects for count outcomes

## Installation
`CountEffects` is currently not on CRAN. The development version of `CountEffects` can be installed directly from this GitHub repository using the additional package `devtools`. Under Windows, please make sure Rtools (http://cran.r-project.org/bin/windows/Rtools/) are installed and no older version of `CountEffects` is currently loaded): 

```
install.packages("devtools")
library(devtools)

install_github("chkiefer/CountEffects")
```

## Run CountEffects

The main function of the package is `countEffects()`. Type `example(countEffects)` for some examples. 
