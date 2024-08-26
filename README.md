#Install


library(devtools)


install_github("daidai905/HCMarker")


#Precautions


#Please check the Idents of your seurat object



HCMarker(c("/public/home/luyq/keti/find_marker/test/test1.rds","/public/home/luyq/keti/find_marker/test/test2.rds"),"ara","root","0.8","0.5",c("e","f"))
