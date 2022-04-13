# When you read in the file automate that the code reads the number of rows then divides the number of rows in 2 so that we can create the vector
library(plyr)
library(crossmatch)

dat <-  read.csv("symmetricMatrix.csv", row.names = 1)
dat
head(dat)
dat <- data.matrix(dat) # convert dataframe to matrix
#isSymmetric.matrix(dat) # check if the matrix is symmetric 

#dat <- Matrix::forceSymmetric(dat) # force the matrix to be symmetric
#head(dat) # show the data
#isSymmetric.matrix(dat)
count(dat,'Tree1')
z <- c(rep(0,20),rep(1,20)) # creates vector
z
crossmatch::crossmatchtest(z,dat) # package :: method
