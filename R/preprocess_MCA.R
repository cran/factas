

preprocess_MCA<-function(df){

n_col<-ncol(df)
hp<-rep(0,nrow(df))



noms<-list()
for (i in 1:n_col){noms[[i]]<-list()}
j<-rep(0,n_col)

for (i in 1:nrow(df))

  {
   num<-0
   k<-1
   while((num==0)&&(k<=j[1]))
      {
      if  (df[i,1]==noms[[1]][[k]]) 
          {num<-k} else {k<-k+1}
      }

if (num==0) {                                 
             j[1]<-j[1]+1                           
             noms[[1]][[j[1]]]<-df[i,1]           
             p<-ncol(df)                      
             df<-data.frame(df[,1:j[1]],hp,df[,(j[1]+1):p]) 
             df[i,1+j[1]]<-1                     
             


} else { df[i,1+k]<-1 } 
                        

for (m in 2:(n_col-1))

{

num<-0
   k<-1
   while((num==0)&&(k<=j[m]))
      {
      if  (df[i,apply(as.matrix((j+1)[1:(m-1)]),2,sum)+1]==noms[[m]][[k]])  
          {num<-k} else {k<-k+1} 
      }



if (num==0) {                               
             j[m]<-j[m]+1                          
             noms[[m]][[j[m]]]<-df[i,apply(as.matrix((j+1)[1:(m-1)]),2,sum)+1]             
              p<-ncol(df)                      
             df<-data.frame(df[,1:(apply(as.matrix((j+1)[1:m]),2,sum)-1)],hp,df[,apply(as.matrix((j+1)[1:m]),2,sum):p]) 
             df[i,apply(as.matrix((j+1)[1:m]),2,sum)]<-1                  
             


}  else { df[i,(apply(as.matrix((j+1)[1:(m-1)]),2,sum)+1+k)]<-1 } 
                       

}



   num<-0
   k<-1
   while((num==0)&&(k<=j[n_col]))
      {
      if  (df[i,apply(as.matrix((j+1)[1:(n_col-1)]),2,sum)+1]==noms[[n_col]][[k]])
          {num<-k} else {k<-k+1}
      }


if (num==0) {
             j[n_col]<-j[n_col]+1
             noms[[n_col]][[j[n_col]]]<-df[i,apply(as.matrix((j+1)[1:(n_col-1)]),2,sum)+1]
             df<-data.frame(df,hp)
             df[i,apply(as.matrix((j+1)[1:n_col]),2,sum)]<-1
             


} else { df[i,apply(as.matrix((j+1)[1:(n_col-1)]),2,sum)+1+k]<-1 }



}


df[,1]<-NULL
for (m in 1:n_col) {df[,apply(as.matrix(j[1:m]),2,sum)+1]<-NULL}
#df

return(df)
}


