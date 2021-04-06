# OTFeatures - Optimal Transport Features  

This is a repository for computing optimal transport features for [brain
morphometry](https://en.wikipedia.org/wiki/Brain_morphometry).

This package contains only the optimal transport feature extraction step and
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
data(brains)
otf = extract.otf.image.3d(brain1, brain2)
par(mfrow=c(2,2))
image(brain1[,,2])
title("From")
image(brain2[,,2])
title("To")
image(brain1[,,2] - brain2[,,2])
title("Difference Per Voxel")
image(otf$difference.from[,,2])
title("Mass Allocation in From")
```

## References 
The optimal transport based morphometry approach was first described in:

> Gerber S, Niethammer M, Styner M, Aylward S. 
> Exploratory Population Analysis with Unbalanced Optimal Transport. 
> Med Image Comput Comput Assist Interv. 2018  
> [Pubmed Link](https://pubmed.ncbi.nlm.nih.gov/31172134/)

A journal article with improvements and additions to the method is in progress.








