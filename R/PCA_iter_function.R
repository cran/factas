
PCA_iter<-function(data,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)

{

 if (stream==TRUE) {iter_init <- 300}

eps <- 0.9

p <- ncol(data)
nb_ind <- nrow(data)


A <- X <- Y <- Corr <- vector("list",nb_fact) 
L1 <- vector(length=nb_fact) 

zbar<-apply(data,2,mean)
C <- ((nb_ind-1)/nb_ind)*cov(data)   
M1 <- 1/diag(C) 

for(i in 1:nb_fact)
   { X[[i]]<-Y[[i]] <- runif(p,min=-1,max=1)
                                    
   }


if (principal_axes==TRUE) {A<-X}
if (eigenvalues==TRUE) 
      {  L1 <- nb_fact:1 
      }
if(corr==TRUE)
      { Corr<-X
      }

n <- 1 

if (stream==TRUE){

while (n<=iter_init) 

   { an <- n^(-eps) 
   
     if (n>=2)
          
          { X <- orth_Gram_Schmidt_metrique_diag(1/M1,Y) 
          }

     for (i in 1:nb_fact)
          
          { temp <- C%*%X[[i]]
            FX <- sum(X[[i]] * temp) / sum(X[[i]]^2 / M1)
            Y[[i]] <- X[[i]] + an*(M1*temp - FX*X[[i]]) 

            if ( principal_axes==TRUE) { A[[i]]=(1/M1)*Y[[i]]
                                      }            

            if (eigenvalues==TRUE) { L1[i] <- L1[i] - an*(L1[i] - FX)
               
                                  } 
            if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(A[[i]]/sqrt(sum(A[[i]]^2*M1)))*sqrt(M1)
                            }
          }

n <- n+1 

}

} else {


debut_chrono<-proc.time()
while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time )) 

   { an <- n^(-eps) 
     

     if (n>=2)
          
          { X <- orth_Gram_Schmidt_metrique_diag(1/M1,Y) 
          }

     for (i in 1:nb_fact)
          
          { temp <- C%*%X[[i]]
            FX <- sum(X[[i]] * temp) / sum(X[[i]]^2 / M1)
            Y[[i]] <- X[[i]] + an*(M1*temp - FX*X[[i]]) 

            if ( principal_axes==TRUE) { A[[i]]=(1/M1)*Y[[i]]
                                      }            

            if (eigenvalues==TRUE) { L1[i] <- L1[i] - an*(L1[i] - FX)  
                                  
                                   } 

            if (corr==TRUE) { Corr[[i]] <- sqrt(L1[i])*(A[[i]]/sqrt(sum(A[[i]]^2*M1)))*sqrt(M1)
                            }
          }


n <- n+1 

}


}

return(list(X=X,L1=L1,A=A,M1=M1,Corr=Corr,zbar=zbar))

}


