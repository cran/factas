
CDA<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,graphics=FALSE,data_init,exec_time,print_step)

{

data<-preprocess_CDA(data)

p<-dim_y <- ncol(data)
dim_r<-groups[1]
dim_s<-groups[2]

X <- A <- Y <-  list()
for (i in 1:nb_fact)
{X[[i]]<-A[[i]]<-Y[[i]]<-c()} 

L1 <- vector(length=nb_fact)


if (stream==TRUE) {

 init <- CDA_iter(data[1:data_init,],groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)
                                                                                                      
                                                                                                       
      X <- Y <- init$X
     L1 <- init$L1
      #C <- init$Cov_Y_calc
      #M <- init$M_calc
      Nn <- init$N_calc
      #G <- init$K_calc%*%init$A2_calc
      A2n <- init$A2_calc                
      A <- init$A
      Gn <- init$M_calc%*%init$K_calc      
      Bn <- init$B_calc
      Corr <- init$Corr

     rbar <- init$rbar
     sbar <- init$sbar
     G1n<-rbind(Gn,t(sbar)-t(rbar)%*%Gn)

     n <- data_init+1
  
table_coord<-matrix(ncol=nb_fact)                              
table_coord<-NULL                                              
x<-c()

debut_chrono<-proc.time()

pas_print<-print_step

eps<-0.9
eps2<-0.6


while ( (n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ) )

   {
     an <- (1/(n+1000000000))^(eps)
     an_vp<-(1/(n+1000000000000))^(eps)
     an2<-(1/(n+100000000))^(eps2)

     y <- t(t(data[n,]))
     r<-t(t(y[1:dim_r]))
     s<-t(t(y[(dim_r+1):dim_y]))
     k<-0
     l<-1 
     while (k==0){if (s[l]==1){k=l} else {l<-l+1}}
     r1<-rbind(r,1)        
     s1<-rbind(s,1)
     rcentre <- r-t(t(rbar))


        if  (!identical(Y,NULL))
      
        { X <- orth_norm_Gram_Schmidt(Nn,Y) 
        }

     if (n%%pas_print==0){writeLines(c(paste("n=",n),"\n"))}    
 
     if (( principal_factors==TRUE)&&(n%%pas_print==0)) {  for (i in 1:nb_fact) {writeLines(c(paste("factor ",i,"\n"),X[[i]],"\n"))}  }     



    for (i in 1:nb_fact)

    {   
          Bn<-Gn%*%A2n
          alpha<-Bn%*%X[[i]]
          dzeta<-Nn%*%X[[i]]
          FX <- (t(alpha)%*%dzeta)/(t(X[[i]])%*%dzeta)
          Y[[i]] <- X[[i]] + an * (alpha - as.numeric(FX)*X[[i]])

          if ((principal_axes==TRUE)||(corr==TRUE)) { A[[i]]<-Nn%*%X[[i]]
                                                      if (n%%pas_print==0) {writeLines(c(paste("axis ",i,"\n"),as.vector(A[[i]]),"\n"))}
                                                    }

          if ((eigenvalues==TRUE)||(corr==TRUE)){ L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                         
                                  if (n%%pas_print==0) {writeLines(c(paste("eigenvalue ",i,"\n"),L1[i],"\n"))}
                                 }



               if (corr==TRUE) {Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*solve(Nn)))))/sqrt(diag(Nn))
                                 if (n%%pas_print==0) {writeLines(c(paste("correlation coefficient ",i,"\n"),Corr[[i]],"\n"))}  
                               }

               if (graphics==TRUE){ x[i]<-t(rcentre)%*%X[[i]]
                                  }

         
    }


     if (graphics==TRUE) {table_coord<-rbind(table_coord,x)}

     Nn <- (Nn+rbar%*%t(rbar))*(1-1/n)+r%*%t(r)/n                       
     rbar<- rbar-(rbar-r)/n 
     sbar<- sbar-(sbar-s)/n      
     Nn<-Nn - rbar%*%t(rbar)


     A2n[k,]<-A2n[k,]-an2*(A2n[k,]-t(r))

     gamma<-r1%*%t(r1)
     G1n<-G1n-an2*gamma%*%G1n
     G1n[,k]<-G1n[,k]+r1
     Gn<-G1n[1:dim_r,]
     


    
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

resultat<-CDA_iter(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)


     X <- resultat$X
     L1 <- resultat$L1
     A <- resultat$A
     Corr <- resultat$Corr
      rbar<-resultat$rbar
      sbar<-resultat$sbar

if (principal_axes==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("axis ",i,"\n"),as.vector(A[[i]]),"\n"))}  }

if (principal_factors==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("factor  ",i,"\n"),X[[i]],"\n"))}  }
    
if (eigenvalues==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("eigenvalue  ",i,"\n"),L1[i],"\n"))}  }     
                                                      
if (corr==TRUE)  {   for (i in 1:nb_fact){writeLines(c(paste("correlation coefficient  ",i,"\n"),Corr[[i]],"\n"))}  }
                                                                                        

if (graphics==TRUE){

palette2<-c()
palette1<-(palette(gray(seq(0,1,len=255))))

for (i in 1:nrow(data)){ palette2[i]<-palette1[floor(i/(floor(nrow(data)/255)+1))+1]}

for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
x<-t(t(data[,1:dim_r])-rbar)%*%X[[i]]
y<-t(t(data[,1:dim_r])-rbar)%*%X[[j]]

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


