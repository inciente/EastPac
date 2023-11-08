# EastPac
Data access and analysis of Eastern Pacific Ocean thermodynamics
Originally created within UCAR and LDEO servers to operate with files stored within. 

Modules: 
- intercomparison.py
- - The most important feature here is the class intercomparison, and the ABC model_run. This helps distribute tasks across multiple model runs and make comparisons between the results of different configurations or ensemble members.
 
- EP_stats.py
- - Includes a suite of mostly statistical functions to be applied to datasets treated in this project. For the most part, the functions here will act through instances of model_run and/or intercomparison.
 
- data_acces.py
- - Access to observational datasets of SST and other reanalyses. Potential development route is to merge this onto a subclass of model_run called SST_analysis, or so. 
