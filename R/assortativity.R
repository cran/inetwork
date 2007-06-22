`assortativity` <-
function(A) {
 k <- numeric()
 for (i in 1:dim(A)[1]) k[i] <- sum(A[i,])
 
 knn <- numeric()
 for (i in 1:dim(A)[1]) {
     knn[i] <- 0
     for (j in 1:dim(A)[1]) knn[i] <- knn[i] + A[i,j]*k[j]
     if (k[i] > 0) knn[i] <- knn[i]/k[i]
     else knn[i] <- 0
 }
 #Knn(k)
 Knnk <- numeric()
 for (h in 1:max(k)) {
     Knnk[h] <- 0
     j <- 0
     for (i in 1:dim(A)[1]) {
         if (k[i] == h) {
            Knnk[h] <- Knnk[h] + knn[i]
            j <- j + 1
         }
     }
     if (j > 0) Knnk[h] <- Knnk[h]/j
 }
 return(list(knn=knn,Knnk=Knnk))
}

