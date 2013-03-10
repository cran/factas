

factas<-function(type,data,groups=NULL,stream=TRUE,nb_fact,principal_factors=TRUE,principal_axes=FALSE,eigenvalues=FALSE,corr=FALSE,graphics=FALSE,data_init,exec_time,print_step)

{

if (type=="PCA"){     
return(PCA(data,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if (type=="MFA"){   
return(MFA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if (type=="GCCA"){       
return(GCCA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if (type=="CCA"){   
return(CCA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if(type=="CDA"){     
return(CDA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if (type=="CA"){      
return(CA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}

else if (type=="MCA"){    
return(MCA(data,groups,stream,nb_fact,principal_factors,principal_axes,eigenvalues,corr,graphics,data_init,exec_time,print_step))
}


else {print("error")}

}



