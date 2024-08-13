workingdir = "" # set working directory before running the script

options(echo=TRUE) # if you want to see commands in output file

args = commandArgs(trailingOnly=T)
rat = as.numeric(args[1])
rat
ind = as.numeric(args[2])
ind

Rcpp::sourceCpp(paste0(workingdir,"/functions/L12proj.cpp"))
Rcpp::sourceCpp(paste0(workingdir,"/functions/SPG_lasso.cpp"))
Rcpp::sourceCpp(paste0(workingdir,"/functions/SPG_glasso.cpp"))
Rcpp::sourceCpp(paste0(workingdir,"/functions/SPG_gbridge.cpp"))
R.files = list.files(path=paste0(workingdir,"/functions"),pattern="*.R",full.names=T)
sapply(R.files,source,.GlobalEnv)

# wavelet parameters
family = "DaubLeAsymm"
fil.num = 4
level = 3

# data generation, 3 groups
dims=c(64,64)
p = prod(dims)
sizes=c(200,200,200)
total = sum(sizes)
cumsize = c(0,cumsum(sizes))
ngroup = length(sizes)

# generate coefficient images
u = seq(0,1,length.out=(dims[1]+1))[-1]
v = seq(0,1,length.out=(dims[2]+1))[-1]
beta = list()
beta[[1]] = matrix(rep(0,prod(dims)),dims[1],dims[2])
beta[[2]] = matrix(rep(0,prod(dims)),dims[1],dims[2])
beta[[3]] = matrix(rep(0,prod(dims)),dims[1],dims[2])
for(a in 1:dims[1]){
  for(b in 1:dims[2]){
    if(u[a]>=1/10 & u[a]<=5/10 & (u[a]/2+v[b])>=7/20 & (v[b]-u[a]/2)<=1/4){beta[[1]][a,b] = 3/5}
    if(u[a]>=5/10 & u[a]<=9/10 & (u[a]/2+v[b])>=19/20 & (v[b]-u[a]/2)<=9/20){beta[[1]][a,b] = 1}
    if(u[a]>=5/40 & u[a]<=19/40 & (u[a]/2+v[b])>=29/80 & (v[b]-u[a]/2)<=19/80){beta[[2]][a,b] = 1}
    if(u[a]>=21/40 & u[a]<=35/40 & (u[a]/2+v[b])>=77/80 & (v[b]-u[a]/2)<=35/80){beta[[2]][a,b] = 1/2}
    if(u[a]>=3/20 & u[a]<=9/20 & (u[a]/2+v[b])>=15/40 & (v[b]-u[a]/2)<=9/40){beta[[3]][a,b] = 4/3}
    if(u[a]>=11/20 & u[a]<=17/20 & (u[a]/2+v[b])>=39/40 & (v[b]-u[a]/2)<=17/40){beta[[3]][a,b] = 4/5}
  }
}

mask = which(c(decomp.2D(list(beta[[1]]),min.level=level))!=0)

X = list()
W = list()
for(m in 1:ngroup){
  for(i in (cumsize[m]+1):cumsize[(m+1)]){
    if(m==1){
      a = rnorm(prod(dims),0,1)
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      a.star = a + noise
      X[[i]] = reconstr.2D(a,dims=dims,min.level=level,family=family,fil.num=fil.num)
      W[[i]] = reconstr.2D(a.star,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }else{
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      X[[i]] = X[[(i-cumsize[m])]]
      W[[i]] = X[[i]] + reconstr.2D(noise,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }
  }
}

g1 = sapply(X[(cumsize[1]+1):cumsize[2]],function(x) sum(c(x*beta[[1]])))
g2 = sapply(X[(cumsize[2]+1):cumsize[3]],function(x) sum(c(x*beta[[2]])))
g3 = sapply(X[(cumsize[3]+1):cumsize[4]],function(x) sum(c(x*beta[[3]])))
g = c(g1,g2,g3)
sigma2 = var(g)/9
Y = g + rnorm(total,0,sqrt(sigma2))

########################
# generate testing data
########################
X.test = list()
W.test = list()
for(m in 1:ngroup){
  for(i in (cumsize[m]+1):cumsize[(m+1)]){
    if(m==1){
      a = rnorm(prod(dims),0,1)
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      a.star = a + noise
      X.test[[i]] = reconstr.2D(a,dims=dims,min.level=level,family=family,fil.num=fil.num)
      W.test[[i]] = reconstr.2D(a.star,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }else{
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      X.test[[i]] = X.test[[(i-cumsize[m])]]
      W.test[[i]] = X.test[[i]] + reconstr.2D(noise,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }
  }
}

g1.test = sapply(X.test[(cumsize[1]+1):cumsize[2]],function(x) sum(c(x*beta[[1]])))
g2.test = sapply(X.test[(cumsize[2]+1):cumsize[3]],function(x) sum(c(x*beta[[2]])))
g3.test = sapply(X.test[(cumsize[3]+1):cumsize[4]],function(x) sum(c(x*beta[[3]])))
g.test = c(g1.test,g2.test,g3.test)
# sigma2 = var(g)/9
Y.test = g.test + rnorm(total,0,sqrt(sigma2)) # use the same sigma2 as training data
var.y = rep(NA,ngroup)
for(k in 1:ngroup){
  var.y[k] = var(Y.test[(cumsize[k]+1):cumsize[(k+1)]])
}




###########################
# generate validation data
###########################
X.vali = list()
W.vali = list()
for(m in 1:ngroup){
  for(i in (cumsize[m]+1):cumsize[(m+1)]){
    if(m==1){
      a = rnorm(prod(dims),0,1)
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      a.star = a + noise
      X.vali[[i]] = reconstr.2D(a,dims=dims,min.level=level,family=family,fil.num=fil.num)
      W.vali[[i]] = reconstr.2D(a.star,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }else{
      noise = rep(0,p)
      if(rat!=0){
        noise[mask] = rnorm(length(mask), 0, 1/sqrt(rat))
      }
      X.vali[[i]] = X.vali[[(i-cumsize[m])]]
      W.vali[[i]] = X.vali[[i]] + reconstr.2D(noise,dims=dims,min.level=level,family=family,fil.num=fil.num)
    }
  }
}


###########
# WNET
###########
library(refund.wave)
library(abind)
beta.wnet = list()
pred.wnet = list()
pmse.wnet = rep(NA,ngroup)
for(k in 1:ngroup){
  fit.wnet = wnet(Y[(cumsize[k]+1):cumsize[(k+1)]],abind(W[(cumsize[k]+1):cumsize[(k+1)]],along=0),
                  min.scale=level,wavelet.family=family,filter.number=fil.num,alpha=c(0.1,0.4,0.7,1))
  print(fit.wnet$tuning.param)
  beta.wnet[[k]] = fit.wnet$fhat
  pred.wnet[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],
                          function(x) sum(c(x*beta.wnet[[k]]))+fit.wnet$coef.params[1])
  pmse.wnet[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.wnet[[k]])^2)/var.y[k]
}
pmse.wnet



###########
# WPCR
###########
beta.wpcr = list()
pred.wpcr = list()
pmse.wpcr = rep(NA,ngroup)
for(k in 1:ngroup){
  fit.wpcr = wcr(Y[(cumsize[k]+1):cumsize[(k+1)]],abind(W[(cumsize[k]+1):cumsize[(k+1)]],along=0),
                 min.scale=level,filter.number=fil.num,wavelet.family=family,
                 nfeatures=200,ncomp=1:20,method='pcr')
  beta.wpcr[[k]] = fit.wpcr$fhat
  pred.wpcr[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],
                          function(x) sum(c(x*beta.wpcr[[k]]))+fit.wpcr$param.coef[1])
  pmse.wpcr[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.wpcr[[k]])^2)/var.y[k]
}
pmse.wpcr







###############################
# glasso w/o noise correction
###############################
library(grpreg)
library(Matrix)
C1 = decomp.2D(W[(cumsize[1]+1):cumsize[2]],min.level=level,family=family,fil.num=fil.num)
C2 = decomp.2D(W[(cumsize[2]+1):cumsize[3]],min.level=level,family=family,fil.num=fil.num)
C3 = decomp.2D(W[(cumsize[3]+1):cumsize[4]],min.level=level,family=family,fil.num=fil.num)
C = as.matrix(bdiag(C1,C2,C3))
cv.fit.grp = cv.grpreg(C,Y,group=rep(1:prod(dims),3),penalty='grLasso',nfolds=5)
fit.grp = grpreg(C,Y,group=rep(1:prod(dims),3),penalty='grLasso',lambda=cv.fit.grp$lambda.min)
beta.glasso = list()
pred.glasso = list()
pmse.glasso = rep(NA,ngroup)
for(k in 1:ngroup){
  beta.glasso[[k]] = reconstr.2D(fit.grp$beta[(2+(k-1)*prod(dims)):(k*prod(dims)+1)],dims=dims,min.level=level,family=family,fil.num=fil.num)
  pred.glasso[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],function(x) sum(c(x*beta.glasso[[k]]))+fit.grp$beta[1])
  pmse.glasso[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.glasso[[k]])^2)/var.y[k]
}

pmse.glasso

lam0 = cv.fit.grp$lambda.min
lam0


######################################
# la4 w/o noise correction
######################################
# several rounds of CV to find optimal params
cv.fit1 = cv.multiSOS(W,Y,family=family,fil.num=fil.num,sizes=sizes,
                      levels=level,lams=lam0*c(0.001,0.01,0.1,1))
cv.fit1$best.lam
cv.fit2 = cv.multiSOS(W,Y,family=family,fil.num=fil.num,sizes=sizes,
                      levels=level,lams=cv.fit1$best.lam*c(0.2,0.5,1,2,5))
cv.fit2$best.lam
cv.fit3 = cv.multiSOS(W,Y,family=family,fil.num=fil.num,sizes=sizes,
                      levels=level,lams=cv.fit2$best.lam*c(1/3,1/2,1,2,3))
cv.fit3$best.lam
fit.gbridge = multiSOS(W,Y,family=family,fil.num=fil.num,
                       min.level=level,sizes=sizes,
                       lam=cv.fit3$best.lam)

beta.gbridge = list()
pred.gbridge = list()
pmse.gbridge = rep(NA,ngroup)
for(k in 1:ngroup){
  beta.gbridge[[k]] = fit.gbridge$beta[[k]]
  pred.gbridge[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],function(x) sum(c(x*beta.gbridge[[k]]))+fit.gbridge$beta0[k])
  pmse.gbridge[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.gbridge[[k]])^2)/var.y[k]
}
pmse.gbridge



####################################
# estimate noise covariance matrix
# with validation samples
####################################
sigma = matrix(0,p,p)
for(i in 1:sizes[1]){
  z1 = c(W.vali[[i]])
  z2 = c(W.vali[[(i + cumsize[2])]])
  z3 = c(W.vali[[(i + cumsize[3])]])

  zbar = (z1+z2+z3)/3
  z = rbind(z1-zbar, z2-zbar, z3-zbar)
  sigma = sigma + (t(z)%*%z)/(2*sizes[1])
}

# get Sigma0 with wavelet transform
B = getB(dims,family,fil.num,level)
Sigma0 = t(B)%*%sigma%*%B

# force Sigma to be diagonal
Sigma = matrix(0,p,p)
diag(Sigma) = diag(Sigma0)

remove(B, sigma, Sigma0)



###################
# noisy lasso
###################
beta.noisy.lasso = list()
pred.noisy.lasso = list()
pmse.noisy.lasso = rep(NA,ngroup)
for(k in 1:ngroup){
  cat("noisy lasso group",k,"\n")
  cv.fit3 = cv.noisy_lasso(X=W[(cumsize[k]+1):cumsize[(k+1)]], Y=Y[(cumsize[k]+1):cumsize[(k+1)]],
                           Sigma=Sigma, zvec=c(1e3), lamvec = lam0*seq(0.1,2,length.out=20),
                           family=family,fil.num=fil.num, min.level=level)
  cat(cv.fit3$best.lam,"\n")
  cat(cv.fit3$best.lam/lam0,"\n","\n")
  fit.noisy.lasso = noisy_lasso(X=W[(cumsize[k]+1):cumsize[(k+1)]], Y=Y[(cumsize[k]+1):cumsize[(k+1)]],
                                Sigma=Sigma, z=1e3, lambda0=cv.fit3$best.lam,
                                family=family,fil.num=fil.num, min.level=level)
  beta.noisy.lasso[[k]] = fit.noisy.lasso$beta
  pred.noisy.lasso[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],
                                 function(x) sum(c(x*beta.noisy.lasso[[k]]))+fit.noisy.lasso$beta0)
  pmse.noisy.lasso[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.noisy.lasso[[k]])^2)/var.y[k]
}
pmse.noisy.lasso


##########################
# noisy glasso
##########################
cv.fit3 = cv.noisy_groupreg(X=W, Y=Y, Sigma=Sigma, sizes=sizes, zvec=c(1e3),
                            lamvec = lam0*seq(0.1,2,length.out=20),
                            method="glasso",
                            family=family,fil.num=fil.num, min.level=level)
cv.fit3$best.lam
cv.fit3$best.lam/lam0

fit.noisy.glasso = noisy_groupreg(X=W, Y=Y, Sigma=Sigma, sizes=sizes, z=1e3,
                                  lambda0=cv.fit3$best.lam,
                                  method="glasso",
                                  family=family,fil.num=fil.num, min.level=level)

beta.noisy.glasso = list()
pred.noisy.glasso = list()
pmse.noisy.glasso = rep(NA,ngroup)
for(k in 1:ngroup){
  beta.noisy.glasso[[k]] = fit.noisy.glasso$beta[[k]]
  pred.noisy.glasso[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],
                                  function(x) sum(c(x*beta.noisy.glasso[[k]]))+fit.noisy.glasso$beta0[k])
  pmse.noisy.glasso[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.noisy.glasso[[k]])^2)/var.y[k]
}
pmse.noisy.glasso




##########################
# noisy gbridge
##########################
cv.fit3 = cv.noisy_groupreg(X=W, Y=Y, Sigma=Sigma, sizes=sizes, zvec=c(1e3),
                            lamvec = lam0*seq(0.1,2,length.out=20),
                            init=fit.noisy.glasso$beta,
                            method="gbridge", maxiter=10, tol0=1e-2,
                            family=family,fil.num=fil.num, min.level=level)
cv.fit3$best.lam
cv.fit3$best.lam/lam0
fit.noisy.gbridge = noisy_groupreg(X=W, Y=Y, Sigma=Sigma, sizes=sizes, z=1e3,
                                   lambda0=cv.fit3$best.lam,
                                   init=fit.noisy.glasso$beta,
                                   method="gbridge", maxiter=20, tol0=1e-2,
                                   family=family,fil.num=fil.num, min.level=level)

beta.noisy.gbridge = list()
pred.noisy.gbridge = list()
pmse.noisy.gbridge = rep(NA,ngroup)
for(k in 1:ngroup){
  beta.noisy.gbridge[[k]] = fit.noisy.gbridge$beta[[k]]
  pred.noisy.gbridge[[k]] = sapply(W.test[(cumsize[k]+1):cumsize[(k+1)]],
                                   function(x) sum(c(x*beta.noisy.gbridge[[k]]))+fit.noisy.gbridge$beta0[k])
  pmse.noisy.gbridge[k] = mean((Y.test[(cumsize[k]+1):cumsize[(k+1)]]-pred.noisy.gbridge[[k]])^2)/var.y[k]
}
pmse.noisy.gbridge


save(beta.wnet,          pmse.wnet,
     beta.wpcr,          pmse.wpcr,
     beta.glasso,        pmse.glasso,
     beta.gbridge,       pmse.gbridge,
     beta.noisy.lasso,   pmse.noisy.lasso,
     beta.noisy.glasso,  pmse.noisy.glasso,
     beta.noisy.gbridge, pmse.noisy.gbridge,
     file=paste0(workingdir,"/results/sim2d_unknown_triangle_error",rat,"_rep",ind,".Rdata"))
