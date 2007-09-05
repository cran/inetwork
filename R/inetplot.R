`inetplot` <-
function(clusters,theta=30,shaft=1,circle=TRUE,singlets=FALSE,labels=TRUE,edges=TRUE,points=TRUE,shadow=TRUE) {
 x <- numeric()
 y <- numeric()
 c <- numeric()
 s <- numeric()
 module <- c(rep(0,length(clusters$singlets)),rep(1:length(clusters$sizes),clusters$sizes))   
 di <- order(clusters$sizes,decreasing=TRUE)
 #theta <- 2*pi/length(di)
 theta <- theta*2*pi/360
 if (theta == 0) theta <- 2*pi/12
 for (i in 1:length(di)) { # no. of communities
     shrink <- clusters$sizes[di[i]]/clusters$sizes[di[1]]
     #shrink <- 1
     centerx <- (i-1)*shaft*shrink*cos((i-1)*theta)
     centery <- (i-1)*shaft*shrink*sin((i-1)*theta)
     #symbols(centerx,centery,circles=shrink,add=TRUE,inches=FALSE)
     tmp1 <- 0
     if (di[i] > 1) tmp1 <- sum(clusters$sizes[1:di[i]-1])
     sizes <- clusters$ringleaders[tmp1 + (1:clusters$sizes[di[i]])]
     size2 <- min(sizes)
     size1 <- max(sizes)    
     tmp <- length(clusters$singlets)
     if (di[i] > 1) tmp <- tmp + sum(clusters$sizes[1:di[i]-1])
     for (j in 1:clusters$sizes[di[i]]) {
         k <- clusters$indices[tmp+j]
         if (circle) {
            phi <- 2*pi/clusters$sizes[di[i]]
            x[k] <- centerx + shrink*cos((j-1)*phi) 
            y[k] <- centery + shrink*sin((j-1)*phi) 
         } else {
            phi <- theta
            x[k] <- centerx + log2(j)*cos((j+2)*phi) 
            y[k] <- centery + log2(j)*sin((j+2)*phi) 
         }
         c[k] <- i+1
         if (size1 != 0) s[k] <- 1.5*clusters$ringleaders[tmp1+j]/size1
         else s[k] <- 1.5*size2
         if (s[k] == 0) s[k] <- 1
     }  
 }
 # for singlets
 if (singlets) {
    radius <- max(max(c(x,y),na.rm=TRUE),abs(min(c(x,y),na.rm=TRUE)))
    radius <- radius + 1
    for (i in 1:length(clusters$singlets)) {
        phi <- 2*pi/length(clusters$singlets)
        k <- clusters$indices[i]
        x[k] <- radius*cos((i-1)*phi) 
        y[k] <- radius*sin((i-1)*phi) 
        c[k] <- 1
        s[k] <- 1
    }  
 }
 op <- par()$mar
 par(mar=c(1.1,1.1,1.1,1.1))
 #plot(x,y,pch=" ",axes=F,ann=F)
 #for (i in 1:length(di)) { # no. of communities
 #    shrink <- 1
 #    centerx <- (i-1)*shaft*shrink*cos((i-1)*theta)
 #    centery <- (i-1)*shaft*shrink*sin((i-1)*theta)
 #    shrink <- clusters$sizes[di[i]]/clusters$sizes[di[1]]
 #    symbols(centerx,centery,circles=shrink,add=T,inches=F)
 #    text(centerx+4,centery,labels=paste("i=",as.character(i)))
 #    #points(centerx,centery,pch=21)
 #}
 plot(x,y,pch=" ",axes=FALSE,ann=FALSE)
 if (edges) {
    for (i in 1:dim(clusters$A)[1]) {
        for (j in i:dim(clusters$A)[1]) {
            if (clusters$A[i,j] == 1) {
               if (module[which(clusters$indices == i)]*module[which(clusters$indices == j)] != 0) {
                  if (module[which(clusters$indices == i)] != module[which(clusters$indices == j)]) { 
                     if (shadow) arrows(x[i],y[i],x[j],y[j],length=0,col="gray")                     else arrows(x[i],y[i],x[j],y[j],length=0)
                  }
               }
            }
        }
    }
 }
 if (edges) {
    for (i in 1:dim(clusters$A)[1]) {
        for (j in i:dim(clusters$A)[1]) {
            if (clusters$A[i,j] == 1) {
               if (module[which(clusters$indices == i)]*module[which(clusters$indices == j)] != 0) {
                  if (module[which(clusters$indices == i)] == module[which(clusters$indices == j)])  
                     arrows(x[i],y[i],x[j],y[j],length=0)
               }
            }
        }
    }
 }

 #symbols(centerx,centery,circles=shrink,add=T,inches=F)
 for (i in 1:dim(clusters$A)[1]) {
     if (labels) {
        if (points) points(x[clusters$indices[i]],y[clusters$indices[i]],pch=19,cex=s[clusters$indices[i]])
        text(x[clusters$indices[i]],y[clusters$indices[i]],labels=c(clusters$singlets,clusters$indivisibles)[i],col=c[clusters$indices[i]],xpd=TRUE)
     } else {
        points(x[clusters$indices[i]],y[clusters$indices[i]],pch=19,col=c[clusters$indices[i]],cex=s[clusters$indices[i]])
     }
 }
 par(mar=op)
}
