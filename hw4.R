
# Load Data
setwd("~/Documents/STAT241/hw4")
data=read.table("./hw4data.data")
head(data)

# Global variables
d=7
N=500
V=1:d
# Graph A
Ga_adj=matrix(c(0,0,0,0,1,1,0,
                0,0,0,0,1,0,0,
                0,0,0,1,0,1,0,
                0,0,0,0,1,0,1,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,1,
                0,0,0,0,0,0,0),d,byrow=T)
Ga_edges=list(c(1,5),c(1,6),c(2,5),c(3,4),c(3,6),c(4,5),c(4,7),
              c(6,7))
# Graph B
Gb_adj=matrix(c(0,1,1,1,0,0,1,
                0,0,0,1,0,1,0,
                0,0,0,1,1,0,1,
                0,0,0,0,0,1,0,
                0,0,0,0,0,0,1,
                0,0,0,0,0,0,0,
                0,0,0,0,0,0,0),d,byrow=T)
Gb_edges=list(c(1,2),c(1,3),c(1,4),c(1,7),c(2,4),c(2,6),c(3,4),c(3,5),
              c(3,7),c(4,6),c(5,7))
# Graph C
Gc_adj=matrix(c(0,1,1,1,1,1,1,
                0,0,1,1,1,1,1,
                0,0,0,1,1,1,1,
                0,0,0,0,1,1,1,
                0,0,0,0,0,1,1,
                0,0,0,0,0,0,1,
                0,0,0,0,0,0,0),d,byrow=T)
Gc_edges=list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6),c(1,7),
              c(2,3),c(2,4),c(2,5),c(2,6),c(2,7),
              c(3,4),c(3,5),c(3,6),c(3,7),c(4,5),c(4,6),
              c(4,7),c(5,6),c(5,7),c(6,7))

# Initialize potentials
init_phi=function(G_adj){
  phi0=array(dim=c(d,d,2,2))
  phi0[,,1,1]=G_adj
  phi0[,,1,2]=G_adj
  phi0[,,2,1]=G_adj
  phi0[,,2,2]=G_adj
  return(phi0)
}


# Joint distribution, non normalized
p0 = function(v,G_edges,phi){
  # v is a d-dimensional vector
  # corresponding to a configuration
  res=1
  # Dual cliques
  for (edge in G_edges){
    i=edge[1]
    j=edge[2]
    res=res*phi[i,j,v[i]+1,v[j]+1]
  }
  return(res)
}


# Possibilities for V
comb = function(n){
  if (n==1)
    return(matrix(c(0,1),2,1))
  else{
    old=comb(n-1)
    add=matrix(NA,2^n,1)
    return(rbind(cbind(matrix(1,2^(n-1),1),old),
                 cbind(matrix(0,2^(n-1),1),old)))
  }
}

# Z
Z = function(G_edges,phi){
  sum(apply(comb(d),1,function(v) p0(v,G_edges,phi)))
}

# Normalized joint distribution
p = function(v,G_edges,phi){
  p0(v,G_edges,phi)/Z(G_edges,phi)
}



# Marginal probability
pC = function(G_edges,phi,C){
  # C=NA,1,NA,0,NA,NA,NA..
  a=which(!is.na(C))[1]
  b=which(!is.na(C))[2]
  m=comb(d)[comb(d)[,a]==C[a], ]
  M=m[m[,b]==C[b], ]
  sum(apply(M,1,function(x) p(x,G_edges,phi)))
}

# Marginal counts
m = function(v,data){
  # v is a d-dimensional vector with NAs
  # indicating which variables to sum
  # in order to have the marginal
  sum(apply(data,1,function(l) all(l==v,na.rm=T)))
}

# IPF
IPF = function(phi0,data,G_edges,eps){
  phi=phi0
  old_phi=phi+2*eps
  it=0
  while (max(abs((old_phi-phi))) > eps){
    old_phi=phi
    # Dual cliques
    for (edge in G_edges){
      a=edge[1]
      b=edge[2]
      for (c in list(c(0,0),c(1,0),c(0,1),c(1,1))){
        v1=c[1]+1
        v2=c[2]+1
        C=rep(NA,d)
        C[a]=c[1]
        C[b]=c[2]
        phi[a,b,v1,v2]=phi[a,b,v1,v2]*m(C,data)/(N*pC(G_edges,phi,C))
      }
    }
    it=+1
    print(max(abs((old_phi-phi))))
  }
  return(list(phi=phi,diff=max(abs((old_phi-phi))),iteration=it))
}
#### Careful here, may want to have a smaller eps #####
phi0a=init_phi(Ga_adj)
IPF_a=IPF(phi0a,data,Ga_edges,0.001)
phi0b=init_phi(Gb_adj)
IPF_b=IPF(phi0b,data,Gb_edges,0.1)
phi0c=init_phi(Gc_adj)
IPF_c=IPF(phi0c,data,Gc_edges,0.1)

IPF_a$phi[3,4,,]
IPF_b$phi[3,4,,]
IPF_c$phi[3,4,,]

# Log likelihood
l = function(G_edges,phi){
  res=0
  for (edge in G_edges){
    a=edge[1]
    b=edge[2]
    for (c in list(c(0,0),c(1,0),c(0,1),c(1,1))){
      v1=c[1]+1
      v2=c[2]+1
      C=rep(NA,d)
      C[a]=c[1]
      C[b]=c[2]
      res=res+m(C,data)*log(phi[a,b,v1,v2])
    }
  }
  res-N*log(Z(G_edges,phi))
}

l_a=l(Ga_edges,IPF_a$phi)
l_b=l(Gb_edges,IPF_b$phi)
l_c=l(Gc_edges,IPF_c$phi)
# Log likelihood of our models
cte=2500
barplot(c(l_a,l_b,l_c)+cte, 
  ylab= paste('Log likelihood +',cte), names.arg = c('A','B','C'),
  main= 'Performance of different graphical models',ylim=c(0,1000))
text(.7,l_a+cte+30, round(l_a))
text(1.9,l_b+cte+30, round(l_b))
text(3.1,l_c+cte+30, round(l_c))


# Mutual Information 
I = function(i,j,data){
  res=0
  for (c in list(c(0,0),c(1,0),c(0,1),c(1,1))){
    # set pij
    v=rep(NA,d)
    v[i]=c[1]
    v[j]=c[2]
    pij=m(v,data)/N
    # set pi
    v[j]=NA
    pi=m(v,data)/N
    # set pj
    v[j]=c[2]
    v[i]=NA
    pj=m(v,data)/N
    res=res+pij*log(pij/(pi*pj))
  }
  return(res)
}


# Set the weights equals to the opposite of I(i,j)
init_weights=function(data){
  m=matrix(NA,d,d)
  for (i in 1:d){
    for (j in 1:d){
      m[i,j]=-I(i,j,data)
    }
  }
  return (m)
}
W=init_weights(data)

# Minimum spanning tree
library(ape)
T_adj=mst(W)
T_adj=T_adj*upper.tri(T_adj)
sum(T_adj*W)

T_edges=list(c(1,2),c(1,7),c(2,6),c(3,7),c(5,7))
phi0T=init_phi(T_adj)
IPF_T=IPF(phi0T,data,T_edges,0.001)
l_T=l(T_edges,IPF_T$phi)

# Summary
barplot(c(l_a,l_b,l_c,l_T)+cte, 
        ylab= paste('Log likelihood +',cte), 
        names.arg = c('A','B','Fully Connected','Max spanning tree'),
        main= 'Performance of different graphical models',ylim=c(0,1000))
text(.7,l_a+cte+30, round(l_a))
text(1.9,l_b+cte+30, round(l_b))
text(3.1,l_c+cte+30, round(l_c))
text(4.3,l_T+cte+30, round(l_T))



