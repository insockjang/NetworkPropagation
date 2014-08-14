networkPropgate_vector<-function(f.0,alpha,Mat){
  # Network progate a signal on a network
  
  # normalized adjacent matrix in column-wise way
  colNorm<-function(mat){
    if(sum(mat)>0){
      return(mat/sum(mat))
    }else{
      return(mat)
    }
  }
  
  AdjMat<-apply(Mat,2,function(x)colNorm(x))
  
  # normalized input seed matrix in colum-wise way for probability vector sum == 1
  
  f.0 <- f.0/sum(f.0)
  
  tof = 1e-12;
  max_iteration = 10000;
  
  if(length(f.0) != nrow(AdjMat)){
    print('Key missmatch error')
    break
  }
  
#   f.t<-f.0
  
  
  itcnt = 1;             
  
  f.prev = f.0;    
  while (itcnt < max_iteration){
    
    f.t = alpha*(AdjMat %*% f.prev) + (1-alpha)*f.0;
    
    
    l1residual = norm(f.prev-f.t ,"1");
    
    f.prev = f.t;
    
    if (l1residual < tof){
      print(paste("Converged at ",itcnt,sep = ""))  
      break
    }
    itcnt = itcnt + 1;        
  }
  
  return(f.t)
}
