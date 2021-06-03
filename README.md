Code for a shiny app for three-level meta-analysis as described in the article 

Assink, M. & Wibbelink, C. J. M. (2016). Fitting three-level meta-analytic models in R: A step-by-step tutorial. The Quantitative Methods for Psychology, 12(3), 154–174. doi:10.20982/tqmp.12.3.p154

The app is free to use and modify under the usual GPL terms. The code is in app.r, for running place the code in a subdirectory threelevel and place a small R script immediately outside the subdirectory with the following two lines of code

library(metafor);library(shiny)

runApp("threelevel")

The R packages metafor, as references in the article cited above, and shiny should be installed in advance. Data should be setup as described in the article, as shown in Table 1, page 157, except that data for each categorical moderator should be included as a single character variable, not as several dummy variables.

March 20., 2021, some minor updates, and included code for funnel plot and modified Egger test as described in the article 

Marengo, D., & Montag, C. (2020). Digital Phenotyping of Big Five Personality via Facebook Data Mining: A Meta-Analysis. Digital Psychology, 1(1), 52–64. https://doi.org/10.24989/dp.v1i1.1823

May 13., 2021, some minor updates, and now possible to get output in Norwegian.

This repo also contains files related to ongoing work together with other researchers at the Centre for Child and Adolescent Mental Health, Eastern and Southern Norway where three-level models are used. At present the file Data for Metaanalysis Emotion Regulation Interventions.csv with data used in a project on emotion regulation is included. When the work is published, more information and a link to the article will be included.
