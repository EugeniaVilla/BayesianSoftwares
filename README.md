
# JAGS STAN and NIMBLE : Hands on comparison
19 February 2020

# Comparison on performance of different Bayesian Softwares: JAGS, STAN and NIMBLE

The repository contains two main elements.

#1 Report

It is available in portable document format(pdf).
The report presents the results of the comparison between three software, Stan, Jags and Nimble, for building analysis methods for statistical models from R, especially for hierarchical models and computationally-intensive methods.
The study is focused on four main model: Linear, Linear-Mixed, Generalized and Acceleration Failure Time.
The main objective of the study is to compare different approaches used by each software and to analyze time and precision of results in order to advice next users that will face this kind of models.
It contains details of Bayesian approach used for each different models analyzed.

#2 MODELS

It is a folder containing a subfolder for each model analyzed, so: LINEAR MODEL, GENERALIZED LINEAR MODEL, LINEAR MIXED MODEL and AFT MODEL subfolders.
Each _NameofModel_ Model folder contains a Final_Code.R, that is the principal file R; different files STAN corresponding to model of different prior's choice of each _NeameofModel_.
In general, the _NameofModel_ Final_Code.R is divided into five main paragraphs:
1. LIBRARIES - PACKAGES
2. DATASET
3. NOT HIERARCHICAL - MODEL 
4. CONJUGATE HIERARCHICAL - MODEL 
5. NOT INFORMATIVE - MODEL 



If you find any error in the code, please let us now (by e-mail addresses reported below).
We will then modify the code, which may be benefical to other readers.


Giulia Gualtieri, giulia.gualtieri@mail.polimi.it
Eugenia Villa, eugenia.villa@mail.polimi.it
Riccardo Vitali, riccardo3.vitali@mail.polimi.it

Department of Matematics,
University of Politecnico di Milano, 
Milano,
Italy

