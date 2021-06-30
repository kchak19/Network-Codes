rm(list=ls())
require(irlba)
#require(igraph)
require(doParallel)
require(Matrix)


n = 20000              #no of vertices
m = 1000             #size of common region
k = 25              #number of subgraphs
b = (n-m)/k
d = 5                #rank of X
repl = 1            #number of replicates

graph = function(P,n)
{
  A = matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] = runif(n*(n-1)/2)
  A = (A + t(A))
  A = (A < P) + 0 ;
  diag(A) = 0
  #A = Matrix(A,sparse=TRUE)
  return(A)
}

ASE = function(A, dim)
{
  A.eig = eigen(A)
  A.eig.values = A.eig$values[1:dim]
  A.eig.vectors = A.eig$vectors[,1:dim]
  if(dim == 1)
    A.coords = sqrt(A.eig.values) * A.eig.vectors
  else
    A.coords = A.eig.vectors %*% diag(sqrt(A.eig.values))
  
  return(A.coords)
}


# ASE = function(A, dim)
# {
#   A.svd = irlba(A, nu = dim, nv = dim)
#   A.svd.values = A.svd$d[1:dim]
#   A.svd.vectors = A.svd$v[,1:dim]
#   if(dim == 1)
#     A.coords = sqrt(A.svd.values) * A.svd.vectors
#   else
#     A.coords = A.svd.vectors %*% diag(sqrt(A.svd.values))
#   
#   return(A.coords)
# }


procr= function(X,Y)
{
  tmp = t(X) %*% Y
  tmp.svd = svd(tmp)
  W = tmp.svd$u %*% t(tmp.svd$v)
  return(list(d = norm(X%*%W - Y, type = "F"), W = W))
}


est = function(A,m,k){          # Function for estimation
  n = nrow(A)
  
  S = apply(A, 1, sum)
  common = tail(order(S),m)
  x = head(order(S),n-m)
  
  samp = c()
  for(i in 1:k)
  {
    samp = rbind(samp, sample(x, size = b, replace = FALSE))
    x = setdiff(x,samp)
  }
  
  ind1 = c(common,samp[1,])
  A.ref = A[ind1,ind1]
  X.hat.ref = ASE(A.ref,d)
  X.ref.0 = X.hat.ref[1:m,]
  X.ref.1 = X.hat.ref[(m+1):(m+b),]
  
  
  
  # This part takes most of the time
  X.part = foreach(i = 2:k, .combine = rbind, .packages = "irlba") %dopar%
    {
      indi = c(common,samp[i,])
      A.sub = A[indi,indi]
      X.hat.sub = ASE(A.sub,d)
      X.sub.0 = X.hat.sub[1:m,]
      
      H = procr(X.sub.0, X.ref.0)$W
      X.hat.sub[(m+1):(m+b),] %*% H
    }
  
  
  
  
  X.fin = matrix(0, n, d)
  
  for(j in 1:length(common)) {X.fin[common[j],] = X.ref.0[j,]}
  for(j in 1:length(samp[1,])) {X.fin[samp[1,j],] = X.ref.1[j,]}
  for(i in 2:k) for(j in 1:length(samp[i,])) {X.fin[samp[i,j],] = X.part[b*(i-2)+j,]}
  
  return(X.fin)
}


registerDoParallel(16)

# N1 = N2 = T1 = T2 = c()


X = matrix(runif(n*d),n,d)
P = X %*% t(X)
P = P/max(P)
diag(P) = 0

A = graph(P,n)


# With subsampling

t3 = proc.time()[3]

X.fin = est(A,m,k)  

P.fin = X.fin %*% t(X.fin)

#P.fin = ifelse(P.fin<0,0,P.fin)
#P.fin = ifelse(P.fin>1,1,P.fin)

t4 = proc.time()[3]

normdiff2 = norm(P-P.fin,type = "F")/norm(P, type= "F")


N2 = normdiff2

T2 = t4-t3



# for(r in 1:repl)
# {
#   X = matrix(runif(n*d),n,d)
#   P = X %*% t(X)
#   P = P/max(P)
#   diag(P) = 0
#   
#   A = graph(P,n)
#   
#   
#   # With subsampling
#   
#   t3 = proc.time()[3]
#   
#   X.fin = est(A,m,k)  
#   
#   P.fin = X.fin %*% t(X.fin)
#   
#   #P.fin = ifelse(P.fin<0,0,P.fin)
#   #P.fin = ifelse(P.fin>1,1,P.fin)
#   
#   t4 = proc.time()[3]
#   
#   normdiff2 = norm(P-P.fin,type = "F")/norm(P, type= "F")
#   
#   N2 = c(N2,normdiff2)
#   T2 = c(T2,t4-t3)
#   
# }

stopImplicitCluster()

# n2 = mean(N2)
# t2 = mean(T2)


result = data.frame(Norm.diff = c(N2), times = c(T2))
rownames(result) = c("with subsamp")


sink("estne2.txt")
print(result)
sink()

