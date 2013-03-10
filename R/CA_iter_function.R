
CA_iter<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
{

if (stream==TRUE) {iter_init <- 1000}

eps <- 0.9

p<-dim_y <- ncol(data)
nb_ind <- nrow(data)
dim_r<-groups[1]
dim_s<-groups[2]


A<-X <- Y <- Corr <- vector("list",nb_fact) 
L1 <- vector(length=nb_fact) 


rbar <- apply(data[,1:dim_r],2,mean)                   
sbar <- apply(data[,((dim_r)+1):dim_y],2,mean)
C <- ((nb_ind-1)/nb_ind)*cov(data)   

Cov_RS<- C[1:dim_r,(dim_r+1):dim_y] 
E_SR <- t(Cov_RS)+sbar%*%t(rbar)
E_RS <- Cov_RS+rbar%*%t(sbar)
Cov_R <- C[1:dim_r,1:dim_r]
Cov_S <- C[(dim_r+1):dim_y,(dim_r+1):dim_y]
E_SS <-Cov_S+sbar%*%t(sbar)
E_RR <-Cov_R+rbar%*%t(rbar)

for(i in 1:nb_fact)
   { X[[i]] <- Y[[i]] <- runif(dim_r,min=-1,max=1)
                                    
   }


N<-diag(E_RR)
M<-1/N


F<-diag(solve(E_SS))*E_SR
F2<-diag(solve(E_RR))*E_RS
G<-E_RS%*%F
B<-M*G


if (principal_axes==TRUE) {A<-X}
if (eigenvalues==TRUE) 
      {  L1 <- nb_fact:1 
      }
if(corr==TRUE) { Corr<-X }
n <- 1 


if (stream==TRUE)

{ 
   while (n<=iter_init) 
   
{ an <- n^(-eps) 

     if (n>=2)
               { X <- orth_Gram_Schmidt_metrique_diag(N,Y)  }

     for (i in 1:nb_fact)

        {
         
          delta<-t(X[[i]])%*%G
          dzeta<-t(X[[i]])*N
          FX <- (delta%*%X[[i]])/(dzeta%*%X[[i]])
          Y[[i]] <- X[[i]] + an * (B%*%X[[i]] - as.numeric(FX)*X[[i]])
         
          if (eigenvalues==TRUE)
                              {L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                     
                              }

          if (principal_axes==TRUE) 
                                { A[[i]]<-N%*%X[[i]] } 

          if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(A[[i]]/sqrt(sum(A[[i]]^2*M)))*sqrt(M)
                            }
        }
    
n <- n+1 

  }

  } else {

debut_chrono<-proc.time()

while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))                

{ an <- n^(-eps)
  

     if (n>=2)          
               { X <- orth_Gram_Schmidt_metrique_diag(N,Y) 
               }         
    for (i in 1:nb_fact)
                       {         
                          gamma<-t(X[[i]])%*%G   
                          dzeta<-t(X[[i]])*N
                          FX <- (gamma%*%X[[i]])/(dzeta%*%X[[i]])
                          Y[[i]] <- X[[i]] + an * (B%*%X[[i]] - as.numeric(FX)*X[[i]])
                          if (eigenvalues==TRUE)
                                             {L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                                                                      
                                             }

                          if (principal_axes==TRUE) 
                                                { A[[i]]<-N*X[[i]] } 
                          
                          if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(A[[i]]/sqrt(sum(A[[i]]^2*M)))*sqrt(M)
                            }

                        }

n <- n+1 ## maj

}
   
}

return(list(N=N,rbar=rbar,sbar=sbar,A=A,F=F,F2=F2,G=G,B=B,E_SR=E_SR,E_SS=E_SS,L1=L1,X=X,Corr=Corr))

}


