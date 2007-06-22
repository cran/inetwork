`icommunity` <-
function(A,nodes=NULL,partite=FALSE) {
 if (dim(A)[1] < 2) return("error: input is not a matrix")
 if (dim(A)[1] != dim(A)[2]) return("error: matrix is not square")
 if (dim(A)[1] < 3) return("network is trivial")
 if (!is.data.frame(nodes) & is.null(nodes)) nodes <- as.character(1:dim(A)[1])

 Q <- 0
 sizes <- numeric()
 indivisibles <- character()
 ringleaders <- numeric()
 arrows <- numeric()
 xcoordinates <- numeric()
 ycoordinates <- numeric()

 B <- matrix(nrow=dim(A)[1],ncol=dim(A)[1],data=0)
 for (i in 1:dim(A)[1]) {
     for (j in 1:dim(A)[1]) {
         B[i,j] = A[i,j]-sum(A[i,])*sum(A[j,])/sum(A)
     }
 }
 if (partite) B <- -B
 eigenB <- eigen(B, symmetric=TRUE)
 iz <- which(eigenB$vectors[,1]==0)
 clusters <- list(Q=Q,sizes=sizes,indivisibles=indivisibles,ringleaders=ringleaders,
arrows=arrows,xcoordinates=xcoordinates,ycoordinates=ycoordinates,indices=iz,singlets=nodes[iz])

 newman <- function(eigenvector1,B,nodes,indices,clusters,origin,singlets) {
  eigenB <- eigen(B, symmetric=TRUE)
  if (eigenB$values[1] <= 0) {
     clusters$sizes <- c(clusters$sizes,dim(B)[1])
     clusters$indivisibles <- c(clusters$indivisibles,nodes)
     clusters$ringleaders <- c(clusters$ringleaders,abs(eigenvector1))
     clusters$indices <- c(clusters$indices,indices)
     for (i in 1:dim(B)[1]) clusters$xcoordinates <- c(clusters$xcoordinates,origin[1])
     for (i in 1:dim(B)[1]) clusters$ycoordinates <- c(clusters$ycoordinates,origin[2]-i)
     return(clusters)
  }
  ip <- which(eigenB$vectors[,1] >0)
  im <- which(eigenB$vectors[,1] <0)
  dQ <- (sum(eigenB$vectors[ip,1]^2) + sum(eigenB$vectors[im,1]^2))*eigenB$values[1]
  clusters$Q <- clusters$Q + dQ
  iz <- which(eigenB$vectors[,1]==0)
  bp <- length(ip)
  bm <- length(im)
  arrows <- origin
  arrows[3] <- origin[1]
  arrows[4] <- origin[2] - dQ
  arrows[5:6] <- arrows[3:4]
  arrows[8] <- arrows[4]
  arrows[9:10] <- arrows[3:4]
  arrows[12] <- arrows[4]
  if (bp >= bm) {
     arrows[7] <- origin[1] + bp
     arrows[11] <- origin[1] - bm    
  } else {
     arrows[7] <- origin[1] + bm 
     arrows[11] <- origin[1] - bp
  }
  
  #print(paste(as.character(dim(B)[1]),as.character(bp),as.character(bm)))
  if (length(iz) != singlets) {
     clusters$sizes <- c(clusters$sizes,rep(1,length(iz)))
     clusters$indivisibles <- c(clusters$indivisibles,nodes[iz])
     clusters$ringleaders <- c(clusters$ringleaders,eigenB$vectors[iz,1])
     clusters$indices <- c(clusters$indices,indices[iz])
     for (i in 1:length(iz)) clusters$xcoordinates <- c(clusters$xcoordinates,arrows[3]+i-1)
     for (i in 1:length(iz)) clusters$ycoordinates <- c(clusters$ycoordinates,arrows[4])
  }

  if (bp == 0 || bm == 0) {
     clusters$sizes <- c(clusters$sizes,bp+bm)
     clusters$indivisibles <- c(clusters$indivisibles,nodes[ip],nodes[im])
     clusters$ringleaders <- c(clusters$ringleaders,abs(eigenB$vectors[ip,1]),abs(eigenB$vectors[im,1]))
     clusters$indices <- c(clusters$indices,indices[ip],indices[im])
     for (i in 1:(bp+bm)) clusters$xcoordinates <- c(clusters$xcoordinates,origin[1])
     for (i in 1:(bp+bm)) clusters$ycoordinates <- c(clusters$ycoordinates,origin[2]-i)
     return(clusters) 
  }

  clusters$arrows <- c(clusters$arrows,arrows)

  if (bp == 1) {
     clusters$sizes <- c(clusters$sizes,1)
     clusters$indivisibles <- c(clusters$indivisibles,nodes[ip])
     clusters$ringleaders <- c(clusters$ringleaders,abs(eigenB$vectors[ip,1]))
     #clusters$Q <- clusters$Q + eigenB$vectors[ip,1]^2 * eigenB$values[1]
     clusters$indices <- c(clusters$indices,indices[ip])
     clusters$xcoordinates <- c(clusters$xcoordinates,arrows[3]-1)
     clusters$ycoordinates <- c(clusters$ycoordinates,arrows[4]-1)
  } else {
     Bp <- matrix(nrow=bp,ncol=bp,data=0)
     for (i in 1:bp) {
         for (j in 1:bp) {
             k <- ip[i]
             l <- ip[j]
             Bp[i,j] <- B[k,l]
         }
     }
     for (i in 1:bp) Bp[i,i] <- Bp[i,i]-sum(Bp[i,])
     #eigenBp <- eigen(Bp, symmetric=TRUE)     
     if (bp >= bm) clusters <- newman(eigenB$vectors[ip,1],Bp,nodes[ip],indices[ip],clusters,arrows[7:8],0)
     else clusters <- newman(eigenB$vectors[ip,1],Bp,nodes[ip],indices[ip],clusters,arrows[11:12],0)
  }

  if (bm == 1) {
     clusters$sizes <- c(clusters$sizes,1)
     clusters$indivisibles <- c(clusters$indivisibles,nodes[im])
     clusters$ringleaders <- c(clusters$ringleaders,abs(eigenB$vectors[im,1]))
     #clusters$Q <- clusters$Q + eigenB$vectors[im,1]^2 * eigenB$values[1]
     clusters$indices <- c(clusters$indices,indices[im])
     clusters$xcoordinates <- c(clusters$xcoordinates,arrows[3]-1)
     clusters$ycoordinates <- c(clusters$ycoordinates,arrows[4]-1)
     return(clusters)
  }
  Bm <- matrix(nrow=bm,ncol=bm,data=0)
  for (i in 1:bm) {
      for (j in 1:bm) {
          k <- im[i]
          l <- im[j]
          Bm[i,j] <- B[k,l]
      }
  }
  for (i in 1:bm)  Bm[i,i] <- Bm[i,i]-sum(Bm[i,])
  #eigenBm <- eigen(Bm, symmetric=TRUE)
  #print(eigenBm$vectors[,1])
  if (bp >= bm) clusters <- newman(eigenB$vectors[im,1],Bm,nodes[im],indices[im],clusters,arrows[11:12],0)
  else clusters <- newman(eigenB$vectors[im,1],Bm,nodes[im],indices[im],clusters,arrows[7:8],0)
 } # end of newman

 if (eigenB$values[1] <= 0) {
    clusters$sizes <- dim(A)[1]-length(iz)
    clusters$indivisibles <- nodes[-iz]
    clusters$ringleaders <- abs(eigenB$vectors[-iz,1])
    clusters$indices <- (1:dim(A)[1])[-iz]
 } else {
    clusters <- newman(eigenB$vectors[,1],B,nodes,1:dim(A)[1],clusters,c(0,0),length(clusters$singlets))
 }
 clusters$Q <- clusters$Q/sum(A)/2
 #if (partite) clusters$Q <- -clusters$Q
 clusters <- c(list(A=A),clusters)
 return(clusters)
} # end of mycluster

