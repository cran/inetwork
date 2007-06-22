`transitivity` <-
function(A) {
 k <- numeric()
 for (i in 1:dim(A)[1]) k[i] <- sum(A[i,])

 c <- numeric()
 for (i in 1:dim(A)[1]) {
     c[i] <- 0
     for (j in 1:dim(A)[1]) {
         for (h in 1:dim(A)[1]) {
             if (h != j && A[i,j] != 0 && A[i,h] != 0 && A[j,h] != 0) {
                c[i] <- c[i] + 1 #A[i,j]*A[i,h]*A[j,h]
             }
         }
     }
     if (k[i] > 1) c[i] <- c[i]/k[i]/(k[i]-1)
     else c[i] <- 0
 }
 #C(k)
 Ck <- numeric()
 for (h in 1:max(k)) {
     Ck[h] <- 0
     j <- 0
     for (i in 1:dim(A)[1]) {
         if (k[i] == h) {
            Ck[h] <- Ck[h] + c[i]
            j <- j + 1
         }
     }
     if (j > 0) Ck[h] <- Ck[h]/j    
 }
 return(list(c=c,Ck=Ck))
}

