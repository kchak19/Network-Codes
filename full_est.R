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



registerDoParallel(16)

# N1 = N2 = T1 = T2 = c()


X = matrix(runif(n*d),n,d)
P = X %*% t(X)
P = P/max(P)
diag(P) = 0

A = graph(P,n)

# Without sub-sampling

t1 = proc.time()[3]

X.hat = ASE(A,d)

P.hat = X.hat %*% t(X.hat)


#P.hat = ifelse(P.hat<0,0,P.hat)
#P.hat = ifelse(P.hat>1,1,P.hat)

t2 = proc.time()[3]

normdiff1 = norm(P-P.hat,type = "F")/norm(P, type="F")


N1 = normdiff1
T1 = t2-t1



# for(r in 1:repl)
# {
#   X = matrix(runif(n*d),n,d)
#   P = X %*% t(X)
#   P = P/max(P)
#   diag(P) = 0
#   
#   A = graph(P,n)
#   
#   # Without sub-sampling
#   
#   t1 = proc.time()[3]
#   
#   X.hat = ASE(A,d)
#   
#   P.hat = X.hat %*% t(X.hat)
#   
#   
#   #P.hat = ifelse(P.hat<0,0,P.hat)
#   #P.hat = ifelse(P.hat>1,1,P.hat)
#   
#   t2 = proc.time()[3]
#   
#   normdiff1 = norm(P-P.hat,type = "F")/norm(P, type="F")
#   
#   
#   N1 = c(N1,normdiff1)
#   T1 = c(T1,t2-t1)
#   
# }

stopImplicitCluster()

# n1 = mean(N1)
# t1 = mean(T1)



result = data.frame(Norm.diff = c(N1), times = c(T1))
rownames(result) = c("without subsamp")


sink("estne1.txt")
print(result)
sink()

