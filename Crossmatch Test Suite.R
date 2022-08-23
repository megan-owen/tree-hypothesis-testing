

 

# Note if you change range for i you have to add more row names
#basic data
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c("a1","Ea1","Va1","dev","pval","approxval")
pathName <-'Desktop/Crossmatch tests/testFolder'
pathName2 <- '/DecimalMatrixGen' # flip the slash if on windows

for (j in 1:9){ # for the folder (may change depending on num folders)
  print(paste0(pathName,j))
  for (i in 1:10){# for the file (may change depending on num files)
    s <- paste0(pathName,j, pathName2,i,'.csv') # format file name string
    #print(s)
    # print(i)
    M <- data.matrix(read.csv(s, row.names = 1)) # turn it to a matrix
    vec <- c(rep(0,(nrow(M)/2)),rep(1,(nrow(M)/2))) # get the vector
    dats <- crossmatchtest(vec,M) #RUN CROSSMATCH TEST
    df[i,]<- c(unlist(dats)) # turn the list into a matrix

  }
  rownames(df) <- c("DecimalMatrixGen1","DecimalMatrixGen2","DecimalMatrixGen3","DecimalMatrixGen4","DecimalMatrixGen5","DecimalMatrixGen6","DecimalMatrixGen7","DecimalMatrixGen8","DecimalMatrixGen9","DecimalMatrixGen10" )
  print(df)

  write.table(df, file = paste0('R_Crossmatch_Results',j,'.csv'), sep = ",",append = FALSE)



}









  