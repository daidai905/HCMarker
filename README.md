**#Install**


library(devtools)


install_github("daidai905/HCMarker")


**#Precautions**


Please check the Idents of your seurat object to ensure that the Idents are clear


**Use displays**


1. If you enter a single seurat object


HCMarker(files ="test1.rds",species ="ara",tissue ="root", percent ="0.8",num="0.5",papers = "PRJNA111")


2.If you enter two or more seurat objects


HCMarker(files=c("test1.rds","test2.rds"),species = "ara",tissue ="root", percent = "0.8", num="0.5",papers = c("PRJNA111","PRJNA222"))



**Input and output interpretation**


**Input: **


files: one or more seurat objects

species: species name


tissue: tissue name 


percent: what percentage of genes you want to keep according to each indicator (pct1, pct2, proportion, avg_Log2FC, COSG_score)


num:  how many genes you want to select with a high confidence score or greater


papers: the name of the subfolder created for each seurat object





