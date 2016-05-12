# Bayesian survival modelling of university outcomes

This repository contains the R code provided as supplementary material for the manuscript "Bayesian survival modelling of university outcomes" which I coauthored with Prof. Mark F. J. Steel (University of Warwick). 

The code implements Bayesian inference for a discrete time competing risks survival model, developed to analyse data related to university outcomes. This model is based on a multinomial logistic formulation. The implementation also includes a Bayesian Model Average approach for variable selection (low dimensional settings only). 

For a full description of the method see
> Bayesian survival modelling of university outcomes. Vallejos, C.A. and Steel, M.F.J. (2016). To appear in *Journal of the Royal Statistical Society - Series A*


## Abstract of the manuscript 

Dropouts and delayed graduations are critical issues in Higher Education Systems worldwide. A key task in this context is to identify risk factors associated to these events, providing potential targets for mitigating policies. For this purpose, we employ a discrete time competing risks survival model, dealing simultaneously with university outcomes and its associated temporal component. We define survival times as the duration of the student's enrolment at university and possible outcomes as graduation or two types of dropout (voluntary and involuntary), exploring the information recorded at admission time (e.g.~educational level of the parents) as potential predictors. While similar strategies have been previously implemented, we extend the previous methods by handling covariate selection within a Bayesian variable selection framework, where model uncertainty is formally addressed through Bayesian Model Averaging. Our methodology is general, however here we focus on undergraduate students enrolled at three selected degree programmes of the Pontificia Universidad Cat\'olica de Chile during the period 2000-2011. Our analysis reveals interesting insights, highlighting the main covariates that influence students' risk of dropout and delayed graduation.
