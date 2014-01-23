
CA<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,graphics=FALSE,data_init,exec_time,print_step)
{

data<-preprocess_CA(data)

p<-dim_y <- ncol(data)
dim_r<-groups[1]
dim_s<-groups[2]

X <- A <- Y <-  list()
for (i in 1:nb_fact)
{X[[i]]<-A[[i]]<-Y[[i]]<-c()} 

L1 <- vector(length=nb_fact)


if (stream==TRUE) {

 init <- CA_iter(data[1:data_init,],groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time) 
                                                                                                                                                                                                             
      X <- Y <- init$X
      L1 <- init$L1
      C <- init$C
      M <- init$M
      Nn <- init$N
      A2n <- init$F                  
      A <- init$A
      A3n<-init$F2
      Bn <- init$B
      Corr<-init$Corr
      rbar <- init$rbar
      sbar <- init$sbar

     n <- data_init+1

table_coord<-matrix(ncol=nb_fact)                       
table_coord<-NULL                                            
x<-c()  

debut_chrono<-proc.time()

pas_print<-print_step

while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))

  {
     an <- (1/(n+100))^(0.9)
     an2 <- (1/(n+100000))^(0.9)
     an_vp <- (1/n)^(0.9)

     y <- t(t(data[n,]))
     r<-t(t(y[1:dim_r]))
     s<-t(t(y[(dim_r+1):dim_y]))
     k_s<-k_r<-0
     l_s<-l_r<-1 
     while (k_s==0){if (s[l_s]==1){k_s=l_s} else {l_s<-l_s+1}}
     while (k_r==0){if (r[l_r]==1){k_r=l_r} else {l_r<-l_r+1}}
     r1<-rbind(r,1)        
     s1<-rbind(s,1)
     rcentre <- r-t(t(rbar))


        if  (!identical(Y,NULL))
      
        { X <- orth_norm_Gram_Schmidt_metrique_diag(Nn,Y)
        }

       if (n%%pas_print==0){writeLines(c(paste("n=",n),"\n"))}    
 
     if (( principal_factors==TRUE)&&(n%%pas_print==0)) {  for (i in 1:nb_fact) {writeLines(c(paste("factor",i,"\n"),X[[i]],"\n"))}  }      


     for (i in 1:nb_fact)

        {
          Bn<-A3n%*%A2n
          alpha<-Bn%*%X[[i]]          
          dzeta<-Nn*X[[i]]
          FX <- (t(alpha)%*%dzeta)/(t(X[[i]])%*%dzeta)
          Y[[i]] <- X[[i]] + an * (alpha - as.numeric(FX)*X[[i]])



          if ((principal_axes==TRUE)||(corr==TRUE)){A[[i]]<-Nn*X[[i]]
                                                    if (n%%pas_print==0) {writeLines(c(paste("axis",i,"\n"),as.vector(A[[i]]),"\n"))}  
                                                   }
          if ((eigenvalues==TRUE)||(corr==TRUE)){L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))        
                                                 if (n%%pas_print==0) {writeLines(c(paste("eigenvalue ",i,"\n"),L1[i],"\n"))}        
                                                 }


               if (corr==TRUE) {Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*(1/Nn)))))/sqrt(Nn)
                                 if (n%%pas_print==0) {writeLines(c(paste("correlation coefficient",i,"\n"),Corr[[i]],"\n"))}  
                               }

               if (graphics==TRUE){ x[i]<-as.numeric(t(rcentre)%*%X[[i]])
                                  }

         
    }


     if (graphics==TRUE) {table_coord<-rbind(table_coord,x)}


rbar<- rbar-(rbar-r)/n 
sbar<- sbar-(sbar-s)/n 

Nn<-(1-1/n)*Nn
Nn[k_r]<-Nn[k_r]+1/n

A2n[k_s,]<-(1-an2)*A2n[k_s,]
A2n[k_s,k_r]<-A2n[k_s,k_r]+an

A3n[k_r,]<-(1-an2)*A3n[k_r,]
A3n[k_r,k_s]<-A3n[k_r,k_s]+an

    
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


 }  else {

resultat<-CA_iter(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,exec_time)


     X <- resultat$X
     L1 <- resultat$L1
     A <- resultat$A
     Corr <- resultat$Corr
      rbar<-resultat$rbar

if (principal_axes==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("axis",i,"\n"),as.vector(A[[i]]),"\n"))}  }

if (principal_factors==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("factor",i,"\n"),X[[i]],"\n"))}  }
    
if (eigenvalues==TRUE) {  for (i in 1:nb_fact){writeLines(c(paste("eigenvalue",i,"\n"),L1[i],"\n"))}  }     
                                                      
if (corr==TRUE)  {   for (i in 1:nb_fact){writeLines(c(paste("correlation coefficient",i,"\n"),Corr[[i]],"\n"))}  }
                                                                                        

if (graphics==TRUE){

palette2<-c()
palette1<-(palette(gray(seq(0,1,len=255))))

for (i in 1:nrow(data)){ palette2[i]<-palette1[floor(i/(floor(nrow(data)/255)+1))+1]}

for (i in 1:(nb_fact-1)){
for (j in (i+1):nb_fact){
x<-t(t(data[1:dim_r])-rbar)%*%X[[i]]
y<-t(t(data[1:dim_r])-rbar)%*%X[[j]]

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









