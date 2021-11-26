Code for a shiny app for three-level meta-analysis as described in the article 

Assink, M. & Wibbelink, C. J. M. (2016). Fitting three-level meta-analytic models in R: A step-by-step tutorial. The Quantitative Methods for Psychology, 12(3), 154–174. doi:10.20982/tqmp.12.3.p154

The app is free to use and modify under the usual GPL terms. It has been thoroughly checked for the main example in the article by Assink and Wibbelink, but it is free software and there is absolutely no warranty. The code is in app2.r or app.r as detailed below, for running place the code in a subdirectory threelevel, rename to app.r if necessary, and place a small R script immediately outside the subdirectory with the following two lines of code

library(metafor);library(shiny)

runApp("threelevel")

The R packages metafor, as referenced in the article cited above, and shiny should be installed in advance. Data should be setup as described in the article, as shown in Table 1, page 157, except that data for each categorical moderator should be included as a single character variable, not as several dummy variables.

March 20., 2021, some minor updates, and included code for funnel plot and modified Egger test as described in the article 

Marengo, D., & Montag, C. (2020). Digital Phenotyping of Big Five Personality via Facebook Data Mining: A Meta-Analysis. Digital Psychology, 1(1), 52–64. https://doi.org/10.24989/dp.v1i1.1823

May 13., 2021, some minor updates, and now possible to get output in Norwegian.

June 7., 2021, detected that changes have been made in the function rma.mv in version 3.0_0 of the R package metafor, that imply that the shiny app needs to be modified to work properly. The modified app is termed app3.r in this repo while the version that works with version 2.4_0 or earlier is termed app2.r. After downloading the preferred version, its name has to be changed to app.r. 

November 5., 2021. Renamed app3.r to app.r . This is the version of the app to be used for rma.mv in version >=3.0_0 of the R package metafor. Only this version of the app is maintained hereafter. In app.r only, not in app2.r, included correlations between overall effect sizes within the categories of a categorical moderator, in analysis with one categorical moderator. Correlations are based on the variance covariance matrix in the rma.mv object for this analysis.

November 26., 2001. This repo also contains files related to ongoing work together with other researchers at the Centre for Child and Adolescent Mental Health, Eastern and Southern Norway where three-level models are used. At present the files R file for Moderator Analysis Emotion Regulation.csv and R file for Overall Synthesis Results Emotion Regulation Metaanalysis.csv with data used in a project on emotion regulation are included. When the work is published, more information and a link to the article will be included.
