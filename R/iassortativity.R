`iassortativity` <-
function(A) {
 if (!is.matrix(A)) return("error: input is not a matrix")
 if (dim(A)[1] != dim(A)[2]) return("error: matrix is not square")
 if (dim(A)[1] < 3) return("network is trivial")
 if (any(A>1)) return("error: not an unweighted adjacency matrix")
 if (any(A<0)) return("error: not an unweighted adjacency matrix")

 ##k <- numeric()
 ##for (i in 1:dim(A)[1]) k[i] <- sum(A[i,])
 k <- rowSums(A)
 
 ##knn <- numeric()
 ##for (i in 1:dim(A)[1]) {
 ##    knn[i] <- 0
 ##    for (j in 1:dim(A)[1]) knn[i] <- knn[i] + A[i,j]*k[j] 
 ##    if (k[i] > 0) knn[i] <- knn[i]/k[i]
 ##    else knn[i] <- 0
 ##}

 f <- function(A,k) return(sum(A*k))
 knn <- apply(A,1,f,k)
 knn <- knn/k
 knn[is.na(knn)] <- 0
  
 #Knn(k)
 ##Knnk <- numeric()
 ##for (h in 1:max(k)) {
 ##    Knnk[h] <- 0
 ##    j <- 0
 ##    for (i in 1:dim(A)[1]) {
 ##        if (k[i] == h) {
 ##           Knnk[h] <- Knnk[h] + knn[i]
 ##           j <- j + 1
 ##        }
 ##    }
 ##    if (j > 0) Knnk[h] <- Knnk[h]/j
 ##}

 ff <- function(i,k,knn) return(sum(knn[which(k == i)])/length(which(k == i)))
 i <- 1:max(k)      
 Knnk <- sapply(i,ff,k,knn)
 Knnk[is.na(Knnk)] <- 0

 return(list(knn=knn,Knnk=Knnk))
}
