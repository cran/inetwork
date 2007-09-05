`ihierarchy` <-
function(clusters,fan=FALSE,spread=1) {
 if (sum(clusters$Q) == 0) return("the network is indivisible")
 op <- par()$mar
 par(mar=c(1.1,2.1,1.1,2.1))

 if (fan) {
    if (spread < 0.01) spread <- 0.01
    x <- abs(clusters$ycoordinates)*cos((spread*clusters$xcoordinates-90)*pi/180)
    y <- abs(clusters$ycoordinates)*sin((spread*clusters$xcoordinates-90)*pi/180)
 } else {
    x <- clusters$xcoordinates
    y <- clusters$ycoordinates
 }

 plot(c(0,x),c(0,y),pch="",axes=FALSE,ann=FALSE)
 for (i in 1:(length(clusters$arrows)/12)) {
     if (fan) {
        x1 <- abs(clusters$arrows[(i-1)*12+2])*cos((spread*clusters$arrows[(i-1)*12+1]-90)*pi/180)
        y1 <- abs(clusters$arrows[(i-1)*12+2])*sin((spread*clusters$arrows[(i-1)*12+1]-90)*pi/180)
        x2 <- abs(clusters$arrows[(i-1)*12+4])*cos((spread*clusters$arrows[(i-1)*12+3]-90)*pi/180)
        y2 <- abs(clusters$arrows[(i-1)*12+4])*sin((spread*clusters$arrows[(i-1)*12+3]-90)*pi/180)
        arrows(x1,y1,x2,y2,length=0)
        r <- abs(clusters$arrows[(i-1)*12+4])
        x1 <- r*cos((spread*clusters$arrows[(i-1)*12+7]-90)*pi/180)
        y1 <- r*sin((spread*clusters$arrows[(i-1)*12+7]-90)*pi/180)
        x2 <- r*cos((spread*clusters$arrows[(i-1)*12+11]-90)*pi/180)
        y2 <- r*sin((spread*clusters$arrows[(i-1)*12+11]-90)*pi/180)
        if (y1 < 0 & y2 < 0) { 
           curve(-sqrt(r*r-x^2),from=-r,to=min(x1,x2),add=TRUE,col="gray",lwd=0.5)
           curve(-sqrt(r*r-x^2),from=max(x1,x2),to=r,add=TRUE,col="gray",lwd=0.5)
           curve(sqrt(r*r-x^2),from=-r,to=r,add=TRUE,col="gray",lwd=0.5)
           curve(-sqrt(r*r-x^2),from=min(x1,x2),to=max(x1,x2),add=TRUE)
        }
        if (y1 >=0 & y2 < 0) { 
           curve(-sqrt(r*r-x^2),from=-r,to=x2,add=TRUE,col="gray",lwd=0.5)
           curve(sqrt(r*r-x^2),from=-r,to=x1,add=TRUE,col="gray",lwd=0.5)
           curve(-sqrt(r*r-x^2),from=x2,to=r,add=TRUE)
           curve(sqrt(r*r-x^2),from=x1,to=r,add=TRUE)  
        }
        if (y1 < 0 & y2 >=0) { 
           curve(-sqrt(r*r-x^2),from=-r,to=x1,add=TRUE,col="gray",lwd=0.5)
           curve(sqrt(r*r-x^2),from=-r,to=x2,add=TRUE,col="gray",lwd=0.5)
           curve(-sqrt(r*r-x^2),from=x1,to=r,add=TRUE)
           curve(sqrt(r*r-x^2),from=x2,to=r,add=TRUE)
        }
        if (y1 >= 0 & y2 >= 0) { 
           curve(sqrt(r*r-x^2),from=-r,to=min(x1,x2),add=TRUE,col="gray",lwd=0.5)
           curve(sqrt(r*r-x^2),from=max(x1,x2),to=r,add=TRUE,col="gray",lwd=0.5)
           curve(-sqrt(r*r-x^2),from=-r,to=r,add=TRUE,col="gray",lwd=0.5)
           curve(sqrt(r*r-x^2),from=min(x1,x2),to=max(x1,x2),add=TRUE)
        }
     } else {
        arrows(clusters$arrows[(i-1)*12+1],clusters$arrows[(i-1)*12+2],clusters$arrows[(i-1)*12+3],
clusters$arrows[(i-1)*12+4],length=0)
        arrows(clusters$arrows[(i-1)*12+7],clusters$arrows[(i-1)*12+8],clusters$arrows[(i-1)*12+11],
clusters$arrows[(i-1)*12+12],length=0)
     }
 }
 nodec <- numeric() 
 for (i in 1:length(clusters$sizes)) nodec <- c(nodec,rep(i,clusters$sizes[i]))
 cols <- order(clusters$sizes,decreasing=TRUE)
 for (i in 1:length(clusters$sizes)) {
     o <- order(clusters$ringleaders[nodec == i],decreasing=TRUE)
     text(x[nodec==i],y[nodec==i],labels=(clusters$indivisibles[nodec==i])[o],col=(which(cols==i)+1),xpd=TRUE)
 }
 par(mar=op)
}

