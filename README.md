Code for a shiny app for three-level meta-analysis as described in the article 

Assink, M. & Wibbelink, C. J. M. (2016). Fitting three-level meta-analytic models in R: A step-by-step tutorial. The Quantitative Methods for Psychology, 12(3), 154–174. doi:10.20982/tqmp.12.3.p154

The code is in app.r, for running place the code in a subdirectory threelevel and place a small R script immediately outside the subdirectory with the following two lines of code

library(metafor);library(shiny)
runApp("threelevel")

The R packages metafor, as references in the article cited above, and shiny should be installed in advance.
