orth_norm_Gram_Schmidt <- function (M,Y)
   
    { 

     nb_fact <- length(Y)
     X <- vector("list",nb_fact)



if (nb_fact==1) { X<-Y} else{

     tempx<-vector("list",nb_fact)
     tempy<-vector("list",nb_fact)
     normX <- vector(length=nb_fact)
     X[[1]] <- as.vector(Y[[1]])
     tempx[[1]]<-t(X[[1]])%*%M
     normX[1] <- as.numeric(tempx[[1]]%*%X[[1]])

     for (i in 2:nb_fact)
         { 


           X[[i]] <- as.vector(Y[[i]])
           tempy[[i]]<-t(as.vector(Y[[i]]))%*%M
           for (j in 1:(i-1))
               { 
           
               X[[i]] <- as.matrix(X[[i]]-(as.numeric((tempy[[i]]%*%X[[j]])/normX[j])) * X[[j]])
               }
          tempx[[i]]<-t(X[[i]])%*%M
          normX[i] <- as.numeric(tempx[[i]]%*%X[[i]])
          }



     for (i in 1:nb_fact)
         { X[[i]] <- X[[i]]/sqrt(normX[i])
         }
  

}
      
     return(X)
     }       

