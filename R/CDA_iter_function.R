
CDA_iter<-function(data,groups,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,exec_time)
{


if (stream==TRUE) {iter_init <- 300}
eps <- 0.9
p<-dim_y <- ncol(data)
nb_ind <- nrow(data)
dim_r<-groups[1]
dim_s<-groups[2]

A <- X <- Y <- Corr <- vector("list",nb_fact) 
L1 <- vector(length=nb_fact)

Cov_Y_calc<- ((nb_ind)/(nb_ind-1))*cov(data)
Cov_RS_calc<-Cov_Y_calc[1:dim_r,(dim_r+1):dim_y]
Cov_SR_calc<-Cov_Y_calc[(dim_r+1):dim_y,1:dim_r]
Cov_R_calc<-Cov_Y_calc[1:dim_r,1:dim_r]
Cov_S_calc<-Cov_Y_calc[(dim_r+1):dim_y,(dim_r+1):dim_y]  

rbar <- apply(data[,1:dim_r],2,mean)                   
sbar <- apply(data[,((dim_r)+1):dim_y],2,mean)

K_calc<-Cov_RS_calc
N_calc<-Cov_R_calc
M_calc<-solve(N_calc)
J_calc<-Cov_SR_calc+sbar%*%t(rbar)
D_calc<- diag(solve(Cov_S_calc + sbar%*%t(sbar)))
A2_calc<-D_calc*J_calc
B_calc<-M_calc%*%K_calc%*%A2_calc



for(i in 1:nb_fact)
   { X[[i]] <- Y[[i]]<- runif(dim_r,min=-1,max=1)
                                   
   }

if (principal_factors==TRUE) {A<-X} 

if (eigenvalues==TRUE) 
      {  L1 <- nb_fact:1 
      }

if(corr==TRUE)
      { Corr<-X
      }

eps<-0.9
n<-1


if (stream==TRUE)

{ 
   while (n<=iter_init)

{  

an <- (100/(n+100))^(eps)


    if  (!identical(Y,NULL))
      
        { X <- orth_Gram_Schmidt(N_calc,Y)
        }


     for (i in 1:nb_fact)

        {
         
          delta<-t(X[[i]])%*%K_calc%*%A2_calc
          dzeta<-t(X[[i]])%*%N_calc
          FX <- (delta%*%X[[i]])/(dzeta%*%X[[i]])
          Y[[i]] <- X[[i]] + an * (B_calc%*%X[[i]] - as.numeric(FX)*X[[i]])
          
          if (eigenvalues==TRUE){ L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                            
                                  }

          if (principal_axes==TRUE)  {A[[i]]<-N_calc%*%X[[i]] } 

          if (corr==TRUE) { Corr[[i]] <- sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M_calc))))/sqrt(diag(N_calc))
                            }
        }
   
n<-n+1

}
} else {

debut_chrono<-proc.time()

while ((n<nrow(data)) && ( (proc.time()-debut_chrono)<exec_time ))             

{

an <- (100/(n+100))^(eps)


    if  (!identical(Y,NULL))
      
        { X <- orth_Gram_Schmidt(N_calc,Y) 
        }

     for (i in 1:nb_fact)

        {
         
          delta<-t(X[[i]])%*%K_calc%*%A2_calc
          dzeta<-t(X[[i]])%*%N_calc
          FX <- (delta%*%X[[i]])/(dzeta%*%X[[i]])
          Y[[i]] <- X[[i]] + an * (B_calc%*%X[[i]] - as.numeric(FX)*X[[i]])

           if (eigenvalues==TRUE){
                                  L1[i] <- L1[i] - an*(L1[i] - as.numeric(FX))
                                         
                                  }

           if (principal_axes==TRUE) { A[[i]]<-N_calc%*%X[[i]] } 

           if (corr==TRUE) { Corr[[i]] <- sqrt(abs(L1[i]))*(as.vector(A[[i]])/sqrt(abs(sum(as.vector(A[[i]])^2*M_calc))))/sqrt(diag(N_calc))
                            }
        }
   
n<-n+1

}


}

return(list(N_calc=N_calc,rbar=rbar,sbar=sbar,A=A,M_calc=M_calc,K_calc=K_calc,J_calc=J_calc,D_calc=D_calc,A2_calc=A2_calc,
B_calc=B_calc,Cov_SR_calc=Cov_SR_calc,Cov_S_calc=Cov_S_calc,L1=L1,X=X,Corr=Corr))

}


