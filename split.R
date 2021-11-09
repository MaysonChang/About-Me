matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}



comMat=function(arr,ori_row,ori_col, r,c)
{
  
  ori=matrix(0,nrow=ori_row,ncol=ori_col)
  
  rg <- (row(ori)-1)%/%r+1
  cg <- (col(ori)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  n=dim(arr)[3]
  for(i in 1:n)
  {
    ori[which(rci==i)]=arr[,,i]
  }
  
  return(ori)
  return(ori)
}