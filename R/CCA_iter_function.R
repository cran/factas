
CCA_iter<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
{

if (stream==TRUE) {iter_init <- 300}

eps <- 0.9

p<-dim_y <- ncol(data)
nb_ind <- nrow(data)
dim_r<-groups[1]
dim_s<-groups[2]

A<-X <- Y <-Corr<- vector("list",nb_fact)
L1 <- vector(length=nb_fact) 

zbar<-apply(data,2,mean)
C <- ((nb_ind-1)/nb_ind)*cov(data)  

Cov_RS<-C[1:dim_r,(dim_r+1):dim_y]
Cov_SR<-C[(dim_r+1):dim_y,1:dim_r]
Cov_R<-C[1:dim_r,1:dim_r]
Cov_S<-C[(dim_r+1):dim_y,(dim_r+1):dim_y]

for(i in 1:nb_fact)
   { X[[i]] <-Y[[i]]<- runif(dim_r,min=-1,max=1)
                                    
   }

rbar <- zbar[1:dim_r]                   
sbar <- zbar[((dim_r)+1):dim_y]

M<-Cov_R
N<-solve(M)

F<-solve(Cov_S)%*%Cov_SR
D<-t(t(rbar))-t(F)%*%t(t(sbar))

G<-Cov_RS%*%F
B<-N%*%G


if (principal_axes==TRUE) {A<-X}
if (eigenvalues==TRUE) 
      {  L1 <- nb_fact:1 
      }
if(corr==TRUE)
      { Corr<-X
      }
n <- 1 

if (stream==TRUE)

{ 
   while (n<=iter_init)
   
{ an <- n^(-eps) 

     if (n>=2)
               { X <- orth_Gram_Schmidt(M,Y)  }

     for (i in 1:nb_fact)

        {          
          gamma<-t(X[[i]])%*%G         
          dzeta<-t(X[[i]])%*%M
          FX <- (gamma%*%X[[i]])/(dzeta%*%X[[i]])
          Y[[i]] <- X[[i]] + an * (B%*%X[[i]] - as.numeric(FX)*X[[i]])
          if (eigenvalues==TRUE)
                              {L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                      
                              }

          if (principal_axes==TRUE) 
                                { A[[i]]<-M%*%X[[i]] 
                                 } 
          if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*N))))/sqrt(diag(M))
                           }
        }
    
n <- n+1 
   }
  } else {

debut_chrono<-proc.time()

while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time )) 
                   

{ an <- n^(-eps) 

     if (n>=2)          
               { X <- orth_Gram_Schmidt(M,Y) 
               }         
    for (i in 1:nb_fact)
                       {         
                          gamma<-t(X[[i]])%*%G   
                          dzeta<-t(X[[i]])%*%M
                          FX <- (gamma%*%X[[i]])/(dzeta%*%X[[i]])
                          Y[[i]] <- X[[i]] + an * (B%*%X[[i]] - as.numeric(FX)*X[[i]])
                          if (eigenvalues==TRUE)
                                             {L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                                                   
                                             }

                          if (principal_axes==TRUE) 
                                                { A[[i]]<-M%*%X[[i]] 
                                                } 
                          if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*N))))/sqrt(diag(M))
                                           }
                        }

n <- n+1 

}
   
}

return(list(A=A,D=D,rbar=rbar,sbar=sbar,F=F,M=M,L1=L1,X=X,N=N,Corr=Corr,zbar=zbar))

}


