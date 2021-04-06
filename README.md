# OTFeatures - Optimal Transport Features  

This is a repository for computing optimal transport features for [brain
morphometry](https://en.wikipedia.org/wiki/Brain_morphometry).

This package contains only te optimal transport feature extraction step and
can be used to integrate optimal transport features with existing morphometry
approaches.

## Install

```R
library(devtools)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
devtools::install_github("samuelgerber/OTFeatures")
```

## Example

```R
library(OTFeatures)
from = array(runif(1000), dim=rep(10, 3))
to = array(runif(1000), dim=rep(10, 3))
otf = extract.otf.image.3d(from, to)
image(otf$difference.from[,,5])
```

## References 
The optimal transport based morphometry approach was first described in:
> Gerber S, Niethammer M, Styner M, Aylward S. 
> Exploratory Population Analysis with Unbalanced Optimal Transport. 
> Med Image Comput Comput Assist Interv. 2018  
> [Pubmed Link](https://pubmed.ncbi.nlm.nih.gov/31172134/)

A journal article with improvements and additions to the method is in progress.








