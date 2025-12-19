# JointMR
JointMR: A joint likelihood-based approach for causal effect estimation in overlapping Mendelian Randomization studies

***Installation***  
`devtools::install_github("hhoulei/JointMR")`  

***Toy Example***  
`library(JointMR)`  
`library(TwoSampleMR)`  
`load("exposure_list_trimmed.Rdata")`  
`load("outcome_list_example.Rdata")`  
`load("outcome_list_all.Rdata")`  
`# The three file can be download from `  
`N_list_example <- list(finn = 486367,`  
`                       MVP  = 432648,`  
`                       EPIC = 22326)`  
`JointMR(exposure_list = exposure_list_trimmed,`  
`        outcome_list = outcome_list_example,`  
`        original_outcome_list = outcome_list_all,`  
`        N_list = N_list_example,`  
`        ancestry = "EUR",`  
`        bootstrap_time = 1000)`  

***Citation***:  
JointMR: A joint likelihood-based approach for causal effect estimation in overlapping Mendelian Randomization studies
Sijia Wu, Lei Hou, Zhongshang Yuan, Xiaoru Sun, Yuanyuan Yu, Hao Chen, Lin Huang,Hongkai Li, Fuzhong Xue

Please contact hou.lei.sph@sdu.edu.cn for any questions. We will continue to update this R package and reduce the problems that may be encountered during its installation.
