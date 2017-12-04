colNorm<-function(mat){
  if(sum(mat)>0){
    return(mat/sum(mat))
  }else{
    return(mat)
  }
}

networkPropgate_vector<-function(f.0,alpha,Mat){
  # Network progate a signal on a network
  
  # normalized adjacent matrix in column-wise way
  
  
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
    # print(itcnt)
  }
  names(f.t)<-names(f.0)
  return(f.t)
}


networkPropgate_randomNetwork<-function(f.0,alpha,Mat,randomIter = 100){
  # Network progate a signal on a random network
  
  f.0 <- f.0/sum(f.0)
  tof = 1e-12;
  max_iteration = 10000;
  
  require(foreach)
  f.tt <- foreach(k = 1:randomIter) %dopar%{
    set.seed(k)
    KK <- sample(1:nrow(Mat),nrow(Mat))
    AdjMat<-apply(Mat[KK, KK],2,function(x)colNorm(x))
    
    # normalized input seed matrix in colum-wise way for probability vector sum == 1
    
    
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
      # print(itcnt)
    }
    
    return(f.t)
    
  }
  # normalized adjacent matrix in column-wise way
  return(f.tt)
}


networkPropgate_randomNetwork1<-function(f.0,alpha,Mat,randomIter = 100){
  # Network progate a signal on a random network
  
  f.0 <- f.0/sum(f.0)
  tof = 1e-12;
  max_iteration = 10000;
  
  f.tt <- vector("list", randomIter)
  
  for(k in 1:randomIter){
    set.seed(k)
    KK <- sample(1:nrow(Mat),nrow(Mat))
    AdjMat<-apply(Mat[KK, KK],2,function(x)colNorm(x))
    
    # normalized input seed matrix in colum-wise way for probability vector sum == 1
    
    
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
      # print(itcnt)
    }
    f.tt[[k]] <- f.t
  }
  
  # normalized adjacent matrix in column-wise way
  return(f.tt)
}


networkPropgate_randomNetwork2 <-function(f.0,alpha,Mat,k){
  # Network progate a signal on a random network
  
  f.0 <- f.0/sum(f.0)
  tof = 1e-12;
  max_iteration = 10000;
  
  
  set.seed(k)
  KK <- sample(1:nrow(Mat),nrow(Mat))
  AdjMat<-apply(Mat[KK, KK],2,function(x)colNorm(x))
  
  # normalized input seed matrix in colum-wise way for probability vector sum == 1
  
  
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
    # print(itcnt)
  }
  return(f.t)
  
}



networkPropgate_randomSeed<-function(f.0,alpha,Mat,randomIter = 100, MC.NUM){
  # Network progate a signal on a random network
  
  f.0 <- f.0/sum(f.0)
  tof = 1e-12;
  max_iteration = 10000;
  AdjMat<-apply(Mat,2,function(x)colNorm(x))
  
  set.seed(k)
  f.0 <- sample(f.0)
  # normalized input seed matrix in colum-wise way for probability vector sum == 1

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
    # print(itcnt)
  }
  return(f.t)
  
}
