library(crossmatch)

dat <-  read.csv("symmetricMatrix.csv", row.names = 1)
# dat # to see the data

dat <- data.matrix(dat) # converts the dataframe into a matrix
 

z <- c(rep(0,(nrow(dat)/2)),rep(1,(nrow(dat)/2))) # creates vector (nrow(dat)/2) 
# gets the number of rows in the matrix then divides that number by 2
# z # to see the vectors
crossmatch::crossmatchtest(z,dat) # package :: method
