
HCMarker <- function(files = files,species =species,tissue = tissue,percent = percent,num = num,papers = papers){
species <<- species
tissue <<- tissue
percent <<- percent
num <<- num
papers <<- papers

  #######################
  high_confidence(files = files,species =species,tissue = tissue,percent = percent,num = num,papers = papers)



sub_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)


if (length(sub_dirs) == 1) {

  sx_file <- list.files(path = sub_dirs[1], pattern = "*sx.txt", full.names = TRUE)
  jj <- read.table(sx_file, sep = '\t', header = TRUE)
} else if (length(sub_dirs) > 1) {
  
  sx_files <- unlist(lapply(sub_dirs, function(dir) list.files(path = dir, pattern = "*sx.txt", full.names = TRUE)))
  jj <- do.call(rbind, lapply(sx_files, function(file) read.table(file, sep = '\t', header = TRUE)))
} else {
  stop("No subdirectories found.")
}

print("CellType signature gene set screening is in progress")





jj <<- unique(jj)





  zj_matrix <- "zj_matrix.txt"
  zj(zj_matrix)


  
  ct_num <- paste0("celltype.num.txt")
  cety_num(zj_matrix, ct_num)


  high <- paste0("score.txt")
  high_score(ct_num, high)


  dayu_0.5 <- paste0(num,".txt")
  unique_gene(high, dayu_0.5)
  file.remove("celltype.num.txt")
}
