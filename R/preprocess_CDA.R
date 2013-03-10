preprocess_CDA<-function(df){

n_col<-ncol(df)
hp<-rep(0,nrow(df))
hp

noms<-list()
j<-0


for (i in 1:nrow(df))

  {
   num<-0
   k<-1
   while((num==0)&&(k<=j))
      {
      if  (df[i,n_col]==noms[[k]])
          {num<-k} else {k<-k+1}
      }

if (num==0) {
             j<-j+1
             noms[[j]]<-df[i,n_col]
             df<-data.frame(df,hp)
             df[i,n_col+j]<-1
             


} else { df[i,n_col+k]<-1 }



}

 
df[,n_col]<-NULL

return(df)

}


