cov2 <- function(x)crossprod(scale(x ,TRUE ,FALSE ))/(NROW(x) -1);cov2 <- compiler::cmpfun(cov2)
library(parallel)
require(doMC)
registerDoMC(cores=7)
library(glmnet)

load('genomatrix.rda')
n = nrow(x)
p = ncol(x)

# non-overlap
b1 = b2 = rep(0,p)
b1[1:50] = rnorm(50)
b2[26:75] = rnorm(50)
# check the true values
maxid = c(1:75)
tam = b1[maxid] %*% ( cov2(x[,maxid]) %*% b2[maxid] )

# simulate the true coefficients
err.mean = err.median = err.la.ols =err.laso =c()
for(ss in 1:50){
  # simulate 2 traits
  y1 = x %*% b1 + rnorm(n)
  y2 = x %*% b2 + rnorm(n)
  
  fit1 = cv.glmnet(x, y1, alpha = 1, parallel = T)
  lambda1 = fit1$lambda.min
  fit2 = cv.glmnet(x, y2, alpha = 1, parallel = T)
  lambda2 = fit2$lambda.min
  
  outlist <- foreach(jj = seq(50) ) %dopar%  {
    foreach(i = seq(2), .packages=c("glmnet")) %dopar%  {
      which <- sample(rep(seq(2), length = n)) == i
      x_y1_1 = x[,predict(glmnet(x[which,], y1[which], alpha= 1, lambda = lambda1), 
                          type = 'nonzero', s = 'lambda.min')[[1]]  ]
      x_y2_1 = x[,predict(glmnet(x[which,], y2[which], alpha= 1, lambda = lambda2), 
                          type = 'nonzero', s = 'lambda.min')[[1]]  ]
      fit1 <- lm(y1[!which] ~ x_y1_1[!which,] + 0 )
      names(fit1$coefficients) <- colnames(x_y1_1)
      bhat1 = coef(fit1)
      bhat1 = na.omit(bhat1)
      fit2 <- lm(y2[!which] ~ x_y2_1[!which,] + 0 )
      names(fit2$coefficients) <- colnames(x_y2_1)
      bhat2 = coef(fit2)
      bhat2 = na.omit(bhat2)
      
      selected.id = unique(c(names(bhat1), names(bhat2) ) ) 
      bh1 = bh2 = matrix(0,nc=1,nr = length(selected.id), dimnames = list(selected.id))
      bh1[names(bhat1),] = bhat1
      bh2[names(bhat2),] = bhat2
      a = NA
      tryCatch({a <-  t(bh1) %*%( cov2(x[,selected.id]) %*%  bh2 )  },  error=function(e){ } )
    }}
  err.mean[ss] = abs(mean(unlist(outlist)) - tam )
  err.median[ss] = abs(median(unlist(outlist)) - tam )
  
  b1.lasso = predict(fit1,type="coef")[-1]
  b2.lasso = predict(fit2,type="coef")[-1]
  sele.lasso.id = unique(c(predict(fit1,type="nonzero")[[1]], predict(fit2,type="nonzero")[[1]] ) ) 
  err.laso[ss] = abs(t(b1.lasso[sele.lasso.id]) %*%( cov2(x[,sele.lasso.id]) %*%  b2.lasso[sele.lasso.id] ) -
                       tam)
  fit.lm.1 <- lm(y1 ~ x[,predict(fit1,type="nonzero")[[1]] ] + 0 )
  names(fit.lm.1$coefficients) <- colnames(x[,predict(fit1,type="nonzero")[[1]] ])
  b1.ols = coef(fit.lm.1)
  b1.ols = na.omit(b1.ols)
  fit.lm.2 <- lm(y2 ~ x[,predict(fit2,type="nonzero")[[1]] ] + 0 )
  names(fit.lm.2$coefficients) <- colnames(x[,predict(fit2,type="nonzero")[[1]] ])
  b2.ols = coef(fit.lm.2)
  b2.ols = na.omit(b2.ols)
  
  sele.id.ols = unique(c(names(b1.ols), names(b2.ols) ) ) 
  bh1 = bh2 = matrix(0,nc=1,nr = length(sele.id.ols), dimnames = list(sele.id.ols))
  bh1[names(b1.ols),] = b1.ols
  bh2[names(b2.ols),] = b2.ols
  
  err.la.ols[ss] = abs(t(bh1) %*%( cov2(x[,sele.id.ols]) %*%  bh2 )-tam )
  print(ss)
}

mean(err.mean);sd(err.mean)
mean(err.median);sd(err.median)
mean(err.laso);sd(err.laso)
mean(err.la.ols);sd(err.la.ols)

quantile(unlist(outlist), prob = c(0.025,0.975))

