norm <- function(A) {
C<-as.matrix(apply(A^(2),2,sum))
D<-(apply(t(C),1,sum))
norme<-sqrt(D)
return(norme)
}











