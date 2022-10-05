# Genetic-Relatedness-HD-Logistic

**Statistical Inference for Genetic Relatedness Based on High-Dimensional Logistic Regression** 

This package contains R functions and simulation codes of a statistical inference method for **genetic relatedness between binary traits based on individual-level genome-wide association data**. Specifically, under the *high-dimensional logistic regression models*, we define parameters characterizing the cross-trait genetic correlation, the genetic covariance and the trait-specific genetic variance. A novel weighted debiasing method is developed for the logistic Lasso estimator  and  computationally efficient debiased estimators are proposed. 


`Simulation-CI.R`, `Simulation-Test.R` and `Simulation-Est.R` contains simple simulation codes for evaluating confidence intervals, statistical tests, and point estimators for genetic correlation.

`main_function` contains the R function of our proposed inference procedure.

Reference: Ma, R., Guo, Z., Cai, T. T., and Li, H. (2022+) Statistical Inference of Genetic Relatedness using High-Dimensional Logistic Regression. Statistica Sinica ([arXiv:2202.10007](https://arxiv.org/pdf/2202.10007.pdf))

For further questions and inquiries, please contact Rong Ma (rongm@stanford.edu).
