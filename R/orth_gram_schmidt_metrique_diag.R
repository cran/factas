
orth_Gram_Schmidt_metrique_diag <- function (M,Y)
   
    { nb_fact <- length(Y)
     X <- vector("list",nb_fact)


if (nb_fact==1) { X<-Y} else{

     normX <- vector(length=nb_fact)
     X[[1]] <- Y[[1]]
     normX[1] <- sum ( X[[1]]^2 * M )

     for (i in 2:nb_fact)
         { X[[i]] <- Y[[i]]
           for (j in 1:(i-1))
               { X[[i]] <- X[[i]]-(sum(Y[[i]] * X[[j]]* M)/normX[j]) * X[[j]]
               }
          normX[i] <- sum (X[[i]]^2 * M)
          }

}        
     return(X)
     }       
