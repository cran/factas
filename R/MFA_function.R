
MFA<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,graphics=FALSE,data_init,exec_time,print_step)

{

p <- ncol(data)
q<-length(groups)
indices_debut_groupes<-c(1,cumsum(groups)+1)


X <- A <- Y <- list()
for (i in 1:nb_fact)
{X[[i]]<-A[[i]]<-Y[[i]]<-c()} 

Xk<-Mk<-Ck<-vector("list",q)
Lk_max<-vector(length=q)

L1 <- vector(length=nb_fact)



if (stream==TRUE) {
 
 init <- MFA_iter(data[1:data_init,],groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)
                                                                                                       
                                                                                                       
      X <- Y <- init$X
     Xk <- init$Xk
     L1 <- init$L1
     Ck <- init$Ck
     Mk <- init$Mk
      A <-init$A
     Corr<-init$Corr

     zbar <- init$zbar 
     
     n <- data_init+1
   


table_coord<-matrix(ncol=nb_fact)                           
table_coord<-NULL                                       
x<-c()

debut_chrono<-proc.time()

pas_print<-print_step


while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))

   {
     an <- 1/n  
     z <- data[n,]
     

   
 M1<-c()
 for (k in 1:q) {
    norme_Xk<-sum(Xk[[k]]^2/Mk[[k]])
    alpha_k <- sum(z[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]*Xk[[k]] )
    beta_k <- sum(zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]*Xk[[k]] ) 
    FXk <- (alpha_k^2 - beta_k^2)/ (norme_Xk)
    Lk_max[k]<- ( t(Xk[[k]])%*%Ck[[k]]%*%Xk[[k]])/ (norme_Xk) 
    Xk[[k]]<-Xk[[k]]+an*(Mk[[k]]*(alpha_k*z[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]-beta_k*zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)])-FXk*Xk[[k]])
    M1<-c(M1,Mk[[k]]/Lk_max[k])
                    }

Var<-1/M1



     if ((!identical(Y,NULL))&(!identical(Var,rep(0,p))))
      
         { X <- orth_norm_Gram_Schmidt_metrique_diag(Var,Y) 
         }

          if (n%%pas_print==0){writeLines(c(paste("n=",n),"\n"))}
         
          if (( principal_factors==TRUE)&&(n%%pas_print==0)) {  for (i in 1:nb_fact) {writeLines(c(paste("factor ",i,"\n"),X[[i]],"\n"))}  }    

  if  (!identical(Var,rep(0,p)))

         { for (i in 1:nb_fact)
        
             { alpha <- sum(z*X[[i]])    
               beta <- sum(zbar*X[[i]]) 
               FX <- (alpha^2-beta^2)/ (sum(X[[i]]^2*Var))
               Y[[i]] <- X[[i]] + an * ((alpha*z-beta*zbar)/Var - FX*X[[i]])
               
               if (( principal_axes==TRUE)||(corr==TRUE)) { A[[i]]=(1/M1)*Y[[i]]
                                                           if (n%%pas_print==0) {writeLines(c(paste("axis ",i,"\n"),A[[i]],"\n"))}  
                                                         }

 
               if ((eigenvalues==TRUE)||(corr==TRUE)) { L1[i] <- L1[i] - an*(L1[i] - FX)
                                                        if (n%%pas_print==0) {writeLines(c(paste("eigenvalue ",i,"\n"),L1[i],"\n"))}     
                                                      }   

               if (corr==TRUE) {Corr[[i]]<-sqrt(abs(L1[i]))*(A[[i]]/sqrt(sum(A[[i]]^2*M1)))*sqrt(M1)
                                 if (n%%pas_print==0) {writeLines(c(paste("correlation coefficient ",i,"\n"),Corr[[i]],"\n"))}  
                               }

               if (graphics==TRUE){ x[i]<-(z-zbar)%*%X[[i]]
                                  }


             }
         }


     if (graphics==TRUE) {table_coord<-rbind(table_coord,x)}

     for (k in 1:q) {
                     Ck[[k]]<- (Ck[[k]]+zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]%*%t(zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]))*
                    (1-1/n)+z[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]%*%t(z[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)])/n
                    }        
      
     zbar <- zbar-(zbar-z)/n             
                     
     for (k in 1:q) {     
                     Ck[[k]]<- Ck[[k]]-zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]%*%t(zbar[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)])
                    Mk[[k]] <-1/diag(Ck[[k]])
                    }

    
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

dev.new()
plot(y,x,xlim=c(-1,1),ylim=c(-1,1),col=palette2,xlab= paste("Axe",i),ylab= paste("Axe",j), main=paste("Representation des invidus dans le plan factoriel",i,j))
axis(2,pos=0)
axis(1,pos=0)
}}

if (corr==TRUE){
for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){

dev.new()
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

} else {


resultat<-MFA_iter(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)


     X <- resultat$X
     L1 <- resultat$L1
     A <- resultat$A
     Corr <- resultat$Corr
     zbar <- resultat$zbar



if (principal_axes==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("axis ",i,"\n"),A[[i]],"\n"))}  }

if (principal_factors==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("factor  ",i,"\n"),X[[i]],"\n"))}  }
    
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

dev.new()
plot(y,x,xlim=c(-1,1),ylim=c(-1,1),col=palette2,xlab= paste("Axe",i),ylab= paste("Axe",j), main=paste("Representation des invidus dans le plan factoriel",i,j))
axis(2,pos=0)
axis(1,pos=0)
}}

if(corr==TRUE){
for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
dev.new()
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


