

#全是基础函数
cosg2 <- function(file = file,file2=file2){
  cosg <- read.table(file,sep = '\t',header = TRUE)


  #提取包含names的列名
  names_cols <- colnames(cosg)[grep("^names", colnames(cosg))]


  names_cols_cleaned <- sub("^names\\.", "", names_cols)



  num_cols <- length(names_cols_cleaned)
  for (i in 1:num_cols){

    col_name <- names_cols_cleaned[i]
    num_rows <- nrow(cosg)



    cosg[[col_name]] <- rep(col_name, num_rows)
  }





  colnames(cosg) <- gsub("\\.", "_", colnames(cosg))



  selected_cols <- colnames(cosg)[!grepl("^names|^scores", colnames(cosg))]


  new_colnames <- sub("^names_|^scores_", "", colnames(cosg))

  colnames(cosg) <- new_colnames



  merged_df <- data.frame()

  for (i in selected_cols){

    df <- cosg[, grepl(paste0("^", i, "\\b"), colnames(cosg))]


    colnames(df) <- c("gene", "cosg", "cell_type_cleaned")





    if (nrow(merged_df) == 0) {
      merged_df <- df
    } else {
      if (ncol(df) != ncol(merged_df)) {
        stop("Number of columns of arguments do not match.")
      }
      merged_df <- rbind(merged_df, df)

    }
  }

  merged_df <- na.omit(merged_df)



  merged_df$cell_type_cleaned <- gsub("\\.", " ", merged_df$cell_type_cleaned)

  write.table(merged_df,file2,sep = '\t',quote =FALSE,row.names =FALSE,col.names=TRUE)
}








score1 <- function(obj=obj,tissue=tissue,article = paper,markertxt =markertxt,filename=filename){
  score_pre <- data.frame(gene= character(), SubSum= numeric(), AllSum=numeric(), Proportion=numeric(), Species=character(), Tissue=character(), Cell_type=character(), DOI=character())

  species <- species

  tissue <- tissue
  article <- article

  obj_markers <- read.table(markertxt,sep = '\t',header = TRUE)
raw_count <- obj[["RNA"]]@data

  for(type in 1:length(unique(obj_markers$cluster))){
    type <- as.character(unique(obj_markers$cluster)[type])
    marker_gene <- obj_markers[obj_markers$cluster == type, "gene"]


    raw_count_subgene <- raw_count[rownames(raw_count) %in% marker_gene, ]

    if (nrow(raw_count_subgene) > 0) {
      raw_count_sum <- as.data.frame(rowSums(raw_count_subgene))
      colnames(raw_count_sum) <- c("AllSum")
      raw_count_sum$gene <- rownames(raw_count_sum)
    } else {
      raw_count_sum <- data.frame(AllSum = numeric(0), gene = character(0))
    }



    sub_result_obj <- subset(obj, idents = type)
    sub_count <- sub_result_obj[["RNA"]]@data
    sub_count <- sub_count[rownames(sub_count) %in% marker_gene, ]
    sub_count_sum <- as.data.frame(rowSums(sub_count))
    colnames(sub_count_sum) <- c("SubSum")
    sub_count_sum$gene <- rownames(sub_count_sum)
    score <- merge(sub_count_sum, raw_count_sum, by = "gene")
    score$Proportion <- score$SubSum / score$AllSum

    score$Species <- rep(species, length(marker_gene))
    score$Tissue <- rep(tissue, length(marker_gene))
    score$Cell_type <- rep(type, length(marker_gene))
    score$Article <- rep(article, length(marker_gene))
    score_pre <- rbind(score_pre, score)
  }

  colnames(score_pre)[1] <- "Marker"
  score <- score_pre[,c("Species", "Tissue", "Cell_type", "Marker",  "Article", "SubSum", "AllSum", "Proportion")]
  colnames(score) <- c("Species", "Tissue", "Cell_type", "Marker",  "Article", "SubSum", "AllSum", "Proportion")
  write.table(score,filename, sep = "\t", row.names = F,quote =FALSE,col.names = TRUE)}





fmfc2 <- function(file=file,scorefile = scorefile,resultfile =resultfile){


  fm <- read.table(file,sep = '\t',header = TRUE)


  colnames(fm) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cell_type","gene")



  sc <- read.table(scorefile,sep = '\t',header = TRUE)
  colnames(sc) <- c("species","tissue","cell_type","gene","article","subsum","allsum","proportion")


  p <- dplyr::inner_join(fm,sc, by = c("gene","cell_type"))

  write.table(p,resultfile,quote = FALSE,sep = '\t',col.names =TRUE,row.names = FALSE)}




ct1_sx <- function(obj = obj,tissue = tissue,paper = paper){


  cg <- read.table("cosg_500.txt",sep = '\t',header =TRUE)

  fm <- read.table("find_marker.txt",sep = '\t',header =TRUE)




  cosg2("cosg_500.txt","recosg.txt")

  co <- read.table("recosg.txt",sep = '\t',header =TRUE)


  score1(obj,tissue,paper,"find_marker.txt","ct1_score.txt")




  fmfc2(file = "find_marker.txt",scorefile = "ct1_score.txt",resultfile ="ct1_score_fm_full.txt")

  full <- read.table("ct1_score_fm_full.txt",sep ='\t',header =TRUE)



  cosg <- read.table("recosg.txt",sep = '\t',header =TRUE)

  colnames(cosg)[3] <- "cell_type"




  all <- dplyr::inner_join(full,cosg,by = c("gene","cell_type"))

  write.table(all,"ct1_jiaoji.txt",quote=FALSE,col.names =TRUE,row.names =FALSE,sep = '\t')



  jj <- read.table("ct1_jiaoji.txt",sep = '\t',header =TRUE)


  jj <- dplyr::arrange(jj,desc(avg_log2FC))

  percent <- as.numeric(percent)
  index <- ceiling(percent * nrow(jj))


  log2 <- jj$avg_log2FC[index]


  jj <- dplyr::arrange(jj,desc(cosg))


  co <-jj$cosg[index]



  jj <- dplyr::arrange(jj,desc(pct.1))

  p1 <- jj$pct.1[index]

  jj <- dplyr::arrange(jj,pct.2)

  p2 <- jj$pct.2[index]


  jj <- dplyr::arrange(jj,desc(proportion))

  index <- ceiling(percent * nrow(jj))

  pro <- jj$proportion[index]


  jj <- jj[jj$avg_log2FC >= log2 &jj$cosg >= co,]



  jj <- jj[jj$pct.1 >= p1 & jj$pct.2 <= p2,]



  jj <- jj[jj$proportion >= pro,]


  f <- paste0(paper,"_sx.txt")
  write.table(jj,f,sep = '\t',col.names =TRUE,row.names =FALSE,quote =FALSE)


file.remove("ct1_score.txt")
file.remove("ct1_score_fm_full.txt")
file.remove("recosg.txt")
file.remove("ct1_jiaoji.txt")
}



zj <- function(filename =filename ){


  result1 <- data.frame(gene = character(),
                        celltype = character(),
                        gene_repeats = numeric(),
                        stringsAsFactors = FALSE
  )

  t <- unique(jj$cell_type)

  for (i in t){
    sub <- jj[jj$cell_type == i,]

    gene_counts <- table(sub$gene)


    sub_result <- data.frame(gene = names(gene_counts),
                             celltype = i,
                             gene_repeats = as.numeric(gene_counts),
                             stringsAsFactors = FALSE)


    result1 <- rbind(result1, sub_result)
  }





 result2 <- reshape2::dcast(result1, gene ~ celltype, value.var = "gene_repeats", fill = 0)


  result2[is.na(result2)] <- 0
  rownames(result2) <- result2$gene
  result2 <- result2[,-1]

  result2$sum <- rowSums(result2)


  result2 <- as.matrix(result2)

  colsum <- colSums(result2)


  result2 <- rbind(result2, sum_col = colsum)


  result2 <- as.data.frame(result2)

  write.table(result2,filename,sep = '\t',col.names =TRUE,row.names =TRUE,quote =FALSE)}




cety_num <- function(zjfilename = zjfilename,cetynum= cetynum){

  zj <- read.table(zjfilename,sep = '\t',header =TRUE)



  zj <- zj[rownames(zj) != "sum_col", ]
  colnames(zj)
  zj <- zj[,colnames(zj) != "sum",]

  zj <- zj[, !colnames(zj) %in% c("Unknown","sum")]

  zj[rownames(zj) == "sum_col",]
  zj[,colnames(zj) == "sum",]




  expre_celltype_num <- numeric(nrow(zj))


  for (i in seq_len(nrow(zj))) {
    gene <- rownames(zj)[i]

    expre_celltype_num[i] <- sum(zj[gene, ] > 0)
  }


  zj$expre_celltype_num <- expre_celltype_num

  write.table(zj,cetynum,sep ='\t',col.names =TRUE,row.names =TRUE,quote =FALSE)}





high_score <- function(zjfilename =zjfilename,scorefile =scorefile){
  zj <- read.table(zjfilename,sep ='\t',header =TRUE)
  res_df <- data.frame(matrix(0, ncol = (ncol(zj) - 1) * 2 + 2, nrow = nrow(zj)))


  colnames(res_df) <- c("Gene", unlist(lapply(unique(colnames(zj)[colnames(zj) != "expre_celltype_num"]),
                                              function(x) c(x, paste0(x, "_score")))), "expre_celltype_num")

  rownames(res_df) <- rownames(zj)

  res_df$Gene <- rownames(zj)

  for (celltype in unique(colnames(zj)[colnames(zj) != "expre_celltype_num"])) {

    a <- zj[, c(celltype, "expre_celltype_num")]


    a <- a[a[, celltype] != 0, ]


    celltype_score <- numeric(nrow(zj))


    for (i in 1:nrow(a)) {
      gene <- rownames(a)[i]
      celltype_score[which(rownames(zj) == gene)] <- (a[i, celltype] - min(a[, celltype]) +1) / (max(a[, celltype]) - min(a[, celltype]) +1) / a[i, "expre_celltype_num"]
    }


    res_df[, celltype] <- zj[, celltype]
    res_df[, paste0(celltype, "_score")] <- celltype_score
  }


  res_df$expre_celltype_num <- zj$expre_celltype_num
  write.table(res_df,scorefile ,sep = '\t',col.names =TRUE,row.names =TRUE,quote =FALSE)

}

unique_gene <- function(scorefile =scorefile,filename =filename){

  results_list <- list()

  res_df <- read.table(scorefile ,sep = '\t',header =TRUE)
  for (b in colnames(res_df[, grepl("score", colnames(res_df))])) {
    print(b)

    num <- as.numeric(num)

    selected_genes <- unique(rownames(res_df[res_df[[b]] >= num, ]))
    print(selected_genes)


    results_list[[b]] <- selected_genes
  }


  max_length <- max(sapply(results_list, length))


  for (b in names(results_list)) {
    length(results_list[[b]]) <- max_length
  }


  results_df <- as.data.frame(results_list)


  results_df[is.na(results_df)] <- ""



  write.table(results_df,filename,sep = '\t',col.names =TRUE,row.names =FALSE,quote =FALSE)}




high_confidence <- function(files = files, species = species, tissue = tissue, percent = percent, num = num, papers = papers) {


  import_seurat_objects <- function(obj_files) {
    print("Starting import_seurat_objects")

    seurat_objects <- list()


    if (is(obj_files, "Seurat")) {
      seurat_objects[[1]] <- obj_files
      print("Single Seurat object provided.")


    } else if (is.character(obj_files)) {
      for (i in seq_along(obj_files)) {
        file_path <- obj_files[i]
        if (!file.exists(file_path)) {
          stop(paste("File does not exist:", file_path))
        }
        print(paste("Reading file:", file_path))
        obj <- tryCatch({
          readRDS(file_path)
        }, error = function(e) {
          stop(paste("Unable to read file:", file_path, "\nError message:", e))
        })
        seurat_objects[[i]] <- obj
      }
      print("Multiple Seurat objects provided.")
    } else {
      stop("Please provide a valid Seurat object or a list of file paths.")
    }

    print("Finished import_seurat_objects")
    return(seurat_objects)
  }

  seurat_objects <- import_seurat_objects(files)


  for (i in seq_along(seurat_objects)) {
    obj <- seurat_objects[[i]]
    paper_name <- papers[i]

    print(paste("Processing object for paper:", paper_name))


    Seurat::DefaultAssay(obj) <- "RNA"


    dir_name <- paste0("./", paper_name)
    dir.create(dir_name, showWarnings = FALSE)
    setwd(dir_name)


    print("Running FindAllMarkers...")
    rds_markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    rds_markers_filtered <- rds_markers[!grepl("Unknown", rds_markers$cluster), ]
    print("Writing find_marker.txt")
    write.table(rds_markers_filtered, "find_marker.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


    print("Running COSG analysis...")
    COSG_markers <- COSG::cosg(obj, groups = 'all', assay = 'RNA', slot = 'data', mu = 1, n_genes_user = 500)
    print("Writing cosg_500.txt")
    write.table(COSG_markers, "cosg_500.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


    print("Calling filter function")
    ct1_sx(obj, tissue, paper = paper_name)


    setwd("..")
  }

  print("All analyses completed.")
  return(invisible(NULL))
}




