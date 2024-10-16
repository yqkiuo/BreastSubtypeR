library(shiny)

options(shiny.maxRequestSize=50*1024^2) 

path="/Users/qiaoyang/OneDrive/OneDrive - Karolinska Institutet/Karolinska Ins/ProjectsAtKI/6.PAM50"
setwd(path)
runApp("ShinyBreastSubtypeR")


## 1. refresh all parameters
## 2. download result format (patientID as the first column)