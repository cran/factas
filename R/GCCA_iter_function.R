
GCCA_iter<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)

{


if (stream==TRUE) {iter_init <- 300}

eps <- 0.9

p <- ncol(data)
nb_ind <- nrow(data)
q<-length(groups)
indices_debut_groupes<-c(1,cumsum(groups)+1)


A<-X <- Y <-Corr<- vector("list",nb_fact) 
M<-Ck<-vector("list",q)
L1 <- vector(length=nb_fact) 

zbar<-apply(data,2,mean)
C <- ((nb_ind-1)/nb_ind)*cov(data)   

for(i in 1:nb_fact)
   { A[[i]] <- Y[[i]] <-runif(p,min=-1,max=1)
                                   
   }


for (k in 1:q){

Ck[[k]] <- C[indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1),indices_debut_groupes[k]:(indices_debut_groupes[k+1]-1)]
M[[k]] <- solve(Ck[[k]])

}

M1<-bdiag(M) 
N1<-solve(M1)
B=C%*%M1

if (principal_factors==TRUE) {X<-A}
if (eigenvalues==TRUE) 
      {  L1 <- nb_fact:1 
      }
if(corr==TRUE)
      { Corr<-A
      }
n <- 1 

if (stream==TRUE)

{ 
   while (n<=iter_init) 

   { 
     an <- n^(-eps) 
   
     if (n>=2)
          
          { A <- orth_Gram_Schmidt(M1,Y) 
          }         
          
     
     for (i in 1:nb_fact)
          
          { 
            temp <- t(A[[i]])%*%M1
            FA <- (temp%*%C%*%t(as.matrix(temp)))/ (temp%*%A[[i]])
            Y[[i]] <- as.matrix( A[[i]] + an * (B%*%A[[i]] - as.numeric(FA)*A[[i]]) )
            

            if ( principal_factors==TRUE) { X[[i]]=N1%*%Y[[i]]
                                      }           

            if (eigenvalues==TRUE) { L1[i] <- L1[i] - an*(L1[i] - as.numeric(FA))  
                                   
                                  } 
            if (corr==TRUE) { Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M1))))/sqrt(abs(diag(N1)))
                            }
          }


n <- n+1 

   }

  } else {

debut_chrono<-proc.time()


while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))            


{ an <- n^(-eps) 
  Y_prec<-as.vector(Y[[nb_fact]])

     if (n>=2)
          
          { A <- orth_Gram_Schmidt(M1,Y) 
          }         
          
     
     for (i in 1:nb_fact)
          
         { 
            temp <- t(A[[i]])%*%M1
            FA <- (temp%*%C%*%t(as.matrix(temp)))/ (temp%*%A[[i]])
            Y[[i]] <- as.matrix(A[[i]] + an * (B%*%A[[i]] - as.numeric(FA)*A[[i]]))
            

            if ( principal_factors==TRUE) { X[[i]]=N1%*%Y[[i]]
                                     }            

            if (eigenvalues==TRUE) { L1[i] <- L1[i] - an*(L1[i] - as.numeric(FA))
                                   
                                  } 
            if (corr==TRUE) { Corr[[i]]<-sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M1))))/sqrt(abs(diag(N1)))
                            }
          }


n <- n+1 

}
   
}

return(list(A=A,C=C,M=M,M1=M1,L1=L1,X=X,Corr=Corr,zbar=zbar))

}




