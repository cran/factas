

preprocess_CA<-function(df){

n_col<-ncol(df)
hp<-rep(0,nrow(df))


noms_r<-noms_s<-list()
j<-0
js<-0


for (i in 1:nrow(df))

  {
   num<-0
   k<-1
   while((num==0)&&(k<=j))
      {
      if  (df[i,1]==noms_r[[k]])
          {num<-k} else {k<-k+1}
      }

if (num==0) {
             j<-j+1
             noms_r[[j]]<-df[i,1]
             p<-ncol(df)
             df<-data.frame(df[,1:j],hp,df[,(j+1):p])
             df[i,1+j]<-1
             


} else { df[i,1+k]<-1 }




   num<-0
   k<-1
   while((num==0)&&(k<=js))
      {
      if  (df[i,j+2]==noms_s[[k]])
          {num<-k} else {k<-k+1}
      }

if (num==0) {
             js<-js+1
             noms_s[[js]]<-df[i,j+2]
             df<-data.frame(df,hp)
             df[i,j+2+js]<-1
             


} else { df[i,j+2+k]<-1 }



}

 
df[,1]<-df[,j+2]<-NULL

return(df)
}


