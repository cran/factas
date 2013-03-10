
GCCA<-function(data,groups,stream=TRUE,nb_fact,principal_factors=FALSE,principal_axes=TRUE,eigenvalues=FALSE,corr=FALSE,graphics=FALSE,data_init,exec_time,print_step)

{

p <- ncol(data)
q<-length(groups)
indices_debut_groupes<-c(1,cumsum(groups)+1)

X <- A <- Y <-  list()
for (i in 1:nb_fact)
{X[[i]]<-A[[i]]<-Y[[i]]<-c()} 

M<-Ck<-vector("list",q)

L1 <- vector(length=nb_fact)


if (stream==TRUE) {
 
 init <- GCCA_iter(data[1:data_init,],groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time) 
                                                                                                      
                                                                                                      
      X <- Y <- init$X
     L1 <- init$L1
     #Ck <- init$Ck
     C <- init$C
     M <- init$M
     M1 <- init$M1
      A <-init$A
      Corr<- init$Corr


     zbar <- init$zbar
     proc_bar<-list(length=q)
     for (i in 1:q){ proc_bar[[i]] <- zbar[indices_debut_groupes[i]:(indices_debut_groupes[i+1]-1)] }

     n <- data_init+1
   

table_coord<-matrix(ncol=nb_fact)                             
table_coord<-NULL                                              
x<-c()

debut_chrono<-proc.time()

pas_print<-print_step

while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))

   {
     an <- 1/n  
     an_vp <- 1/n
     z <- t(t(data[n,]))
     Y_prec<-as.vector(Y[[nb_fact]])
     
   
     if  (!identical(Y,NULL))
      
        { A <- orth_norm_Gram_Schmidt(M1,Y) 
        }

    if (n%%pas_print==0){writeLines(c(paste("n=",n),"\n"))}    
 
     if (( principal_axes==TRUE)&&(n%%pas_print==0)) {  for (i in 1:nb_fact) {writeLines(c(paste("axis ",i,"\n"),A[[i]],"\n"))}  }     
  

    for (i in 1:nb_fact)
        
        { B<-C%*%M1
          
          temp1 <- t(A[[i]])%*%M1
          
          FA <- (temp1%*%C%*%t(temp1)) / (temp1%*%A[[i]])
          Y[[i]] <- A[[i]] + an * (C%*%t(temp1) - as.numeric(FA)*A[[i]])

          if (principal_factors==TRUE) { N1<-solve(M1)
                                         X[[i]]<-N1%*%A[[i]] }
                                         if (n%%pas_print==0) {writeLines(c(paste("factor ",i,"\n"),as.vector(X[[i]]),"\n"))}

          if ((eigenvalues==TRUE)||(corr==TRUE)) {L1[i] <- L1[i] - an_vp*(L1[i] - as.numeric(FA))
                                                           
                                                           if (n%%pas_print==0) {writeLines(c(paste("eigenvalue ",i,"\n"),L1[i],"\n"))}
                                                 }
          
         if (corr==TRUE) { if (principal_factors==TRUE) {Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M1))))/sqrt(abs(diag(N1)))
                                                      } else {Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M1))))/sqrt(abs(diag(solve(M1))))}
                                                       if (n%%pas_print==0) {writeLines(c(paste("correlation coefficient ",i,"\n"),Corr[[i]],"\n"))}  
                         }

               if (graphics==TRUE){ x[i]<-as.numeric(t(z-zbar)%*%X[[i]])
                                  }


         
    }



if (graphics==TRUE) {table_coord<-rbind(table_coord,x)}

C <- (C+zbar%*%t(zbar))*(1-1/n)+z%*%t(z)/n 
zbar<-zbar-(zbar-z)/n                 
C <- C - zbar%*%t(zbar)

for (i in 1:q)   {

proc<-z[indices_debut_groupes[i]:(indices_debut_groupes[i+1]-1)]
proc_centre<-proc-proc_bar[[i]]

gammacentre<-t(proc_centre)%*%M[[i]]                   
phi<-t(gammacentre)%*%gammacentre       
rho<-gammacentre%*%proc_centre
s<-((n+1)/(n*(n+1+as.numeric(rho))))                   
M[[i]] <- ((n+1)/n)*M[[i]] - as.numeric(s)*phi      

proc_bar[[i]]<- proc_bar[[i]]-(proc_bar[[i]]-proc)/n

}
 
M1<-bdiag(M) 

    
n <- n+1 
    
   }

if (graphics==TRUE) {

palette2<-c()
palette1<-(palette(gray(seq(0,.9,len=255))))
for (i in 1:length(table_coord[,1])){ palette2[i]<-palette1[floor(i/(floor(nrow(data)/255)+1))+1]}

for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
x<-table_coord[,i]      
y<-table_coord[,j]     

x_max<-max(max(abs(x)),1)
y_max<-max(max(abs(y)),1)

x<-x/x_max
y<-y/y_max

x11()
plot(y,x,xlim=c(-1,1),ylim=c(-1,1),col=palette2,xlab= paste("Axe",i),ylab= paste("Axe",j), main=paste("Representation des invidus dans le plan factoriel",i,j))
axis(2,pos=0)
axis(1,pos=0)
}}

if (corr==TRUE){
for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){

x11()
par(pty="s")
plot(as.vector(Corr[[i]]),as.vector(Corr[[j]]), xlim=c(-1,1), ylim=c(-1,1),pch="+",xlab=paste("Axe",i),ylab=paste("Axe",j),main=paste("Cercle des correlations des variables avec les facteurs",i,"et",j))
text(as.vector(Corr[[i]])+0.06,as.vector(Corr[[j]]), labels=1:ncol(data))
symbols(0,0,circles=1, inches=FALSE, add=TRUE)
axis(2,pos=0)
axis(1,pos=0)
}
}

}


}



 }  else {

resultat<-GCCA_iter(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)


     X <- resultat$X
     L1 <- resultat$L1
     A <- resultat$A
     Corr <- resultat$Corr
     zbar<-resultat$zbar


if (principal_axes==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("axis ",i,"\n"),A[[i]],"\n"))}  }

if (principal_factors==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("factor  ",i,"\n"),as.vector(X[[i]]),"\n"))}  }
    
if (eigenvalues==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("eigenvalue  ",i,"\n"),L1[i],"\n"))}  }     
                                                      
if (corr==TRUE)  {   for (i in 1:nb_fact){writeLines(c(paste("correlation coefficient  ",i,"\n"),Corr[[i]],"\n"))}  }
                                                                                        

if (graphics==TRUE){

palette2<-c()
palette1<-(palette(gray(seq(0,1,len=255))))

for (i in 1:nrow(data)){ palette2[i]<-palette1[floor(i/(floor(nrow(data)/255)+1))+1]}

for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
x<-t(t(data)-zbar)%*%X[[i]]
y<-t(t(data)-zbar)%*%X[[j]]

x_max<-max(max(abs(x)),1)
y_max<-max(max(abs(y)),1)

x<-x/x_max
y<-y/y_max

x11()
plot(y,x,xlim=c(-1,1),ylim=c(-1,1),col=palette2,xlab= paste("Axe",i),ylab= paste("Axe",j), main=paste("Representation des invidus dans le plan factoriel",i,j))
axis(2,pos=0)
axis(1,pos=0)
}}

if(corr==TRUE){
for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
x11()
par(pty="s")
plot(as.vector(Corr[[i]]),as.vector(Corr[[j]]), xlim=c(-1,1), ylim=c(-1,1),pch="+",xlab=paste("Axe",i),ylab=paste("Axe",j),main=paste("Cercle des correlations des variables avec les facteurs",i,"et",j))
text(as.vector(Corr[[i]])+0.06,as.vector(Corr[[j]]), labels=1:ncol(data))
symbols(0,0,circles=1, inches=FALSE, add=TRUE)
axis(2,pos=0)
axis(1,pos=0)
}}
}

}

}


}



