**#Install**


library(devtools)


install_github("daidai905/HCMarker")


**#Precautions**


Please check the Idents of your seurat object to ensure that the Idents are clear


Use displays


1. If you enter a single seurat object


HCMarker(files ="test1.rds",species ="ara",tissue ="root", percent ="0.8",num="0.5",papers = "PRJNA111")


2.If you enter two or more seurat objects


HCMarker(files=c("test1.rds","test2.rds"),species = "ara",tissue ="root", percent = "0.8", num="0.5",papers = c("PRJNA111","PRJNA222"))



Input vs. output interpretation

