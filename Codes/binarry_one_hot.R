# This function is for converting a fasta file into binary
# one-hot encoded .csv file. 
# Input: 
# 1) data =  file name or location
# 2) format = format of the inputfile (default: fasta, phylip)


binary_one_hot <- function(data, file_format = "fasta", output_file){

  ################# install required packages ##############
  
  if (!requireNamespace("seqinr", quietly = TRUE))
    install.packages("seqinr")
	
  if (!requireNamespace("caret", quietly = TRUE))
    install.packages("caret")

  if (!requireNamespace("stringr", quietly = TRUE))
    install.packages("stringr")	
  
  ################# check required packages ##############
  if (!library('seqinr',logical.return = TRUE)){
    stop("'seqinr' package not found, please install it to run binary_one_hot")
  }
  if (!library('caret',logical.return = TRUE)){
    stop("'caret' package not found, please install it to run binary_one_hot")
  }
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run binary_one_hot")
  }
  
  a <- seqinr::read.alignment(data, format = file_format, forceToLower = FALSE)
  b <- a$seq[12:35]
  b<- lapply(b, function(x) unlist(strsplit(unlist(x), "")) )
  df <- data.frame(matrix(unlist(b), nrow=length(b), byrow=T),stringsAsFactors=FALSE)
  df <- df[vapply(df, function(x) length(unique(x)) > 1, logical(1L))]
  x1 <- caret::dummyVars(~. , data = df)
  df <- data.frame(predict(x1, newdata = df))
  df <- df[,!str_detect(colnames(df), "\\.")]
  rownames(df) <- a$nam[12:35]
  x <- str_remove(colnames(df), "X")
  x1 <- gsub("\\D","",x)
  x2 <- toupper(gsub("[[:digit:]]","",x))
  x  <- paste(x1, x2, sep = "_")
  x2 <- paste("C",str_remove(unlist(str_split(basename(data), "[.]"))[1], "_C12"), sep = "")
  colnames(df) <- paste(x2, x, sep = "_")
  write.csv(df, file = paste(x2, "_C12.csv", sep = ""))
  #rm_row <- apply(df, 2, function(x){length(unique(x))})
  #write.csv(df[,-as.integer(which(rm_row == 1))], file = paste(x2, "_C12.csv", sep = ""))
  rm(df)
  rm(a)
}
