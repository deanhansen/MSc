# Author(s): Dean Hansen

# Posterior_Ratemaking_2017_Table <-
#   matrix(c(24408,1916,296,69,12,6,0,
#            1068,317,61,21,6,2,2,
#            203,71,18,6,2,1,1,
#            49,14,8,3,3,1,0,
#            11,6,2,0,1,0,0,
#            2,0,0,0,0,0,1,
#            1,0,0,1,0,0,0,
#            0,0,0,0,0,0,0,
#            0,0,1,0,0,0,0),
#          nrow=9,ncol=7,byrow=TRUE)

Posterior_Ratemaking_2017_Table <- 
  matrix(c(71087,3722,807,219,51,14,4,0,
           3022,686,184,71,26,10,3,1,
           574,138,55,15,8,4,1,1,
           149,42,21,6,6,1,0,1,
           29,15,3,2,1,1,0,0,
           4,1,0,0,0,0,2,0,
           2,1,0,1,0,0,0,0,
           1,0,0,1,0,0,0,0,
           0,0,1,0,0,0,0,0),
         nrow=9,ncol=8,byrow=TRUE)

# Define dimensions
num_rows <- sum(Posterior_Ratemaking_2017_Table)
num_cols <- 2

# Create an empty matrix to store expanded data
Posterior_Ratemaking_2017_Rows <- matrix(0, nrow = num_rows, ncol = num_cols)

# Populate the expanded matrix
current_row <- 1
for (i in 1:nrow(Posterior_Ratemaking_2017_Table)) {
  for (j in 1:ncol(Posterior_Ratemaking_2017_Table)) {
    count <- Posterior_Ratemaking_2017_Table[i, j]
    if (count > 0) {
      im <- i-1
      jm <- j-1
      Posterior_Ratemaking_2017_Rows[current_row:(current_row + count - 1), 1] <- im
      Posterior_Ratemaking_2017_Rows[current_row:(current_row + count - 1), 2] <- jm
      current_row <- current_row + count
    }
  }
}

write.table(Posterior_Ratemaking_2017_Table, "Main/Data/Bermudez and Karlis_2017/Posterior_Ratemaking_2017_Table.csv",quote=FALSE,row.names=0:8,col.names=0:7,sep=",")
write.table(Posterior_Ratemaking_2017_Rows, "Main/Data/Bermudez and Karlis_2017/Posterior_Ratemaking_2017_Rows.csv",quote=FALSE,row.names=FALSE,col.names=c("N1","N2"),sep=",")

