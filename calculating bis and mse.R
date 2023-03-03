library(parallel)
library(MASS)

confident_interval_fun = function(beta_mat){
  mean_beta = apply(beta_mat,2,mean)
  upper_beta = 0
  lower_beta = 0
  p = dim(beta_mat)[2]
  for(i in c(1:p)){
    upper_beta[i] = quantile(beta_mat[,i],0.975)
    lower_beta[i] = quantile(beta_mat[,i],0.025)
  }
  beta_co = rbind(lower_beta,upper_beta)
}
est_beta_fun = function(n){
  newton_lam = function(c0){
    xt = c0; k0 = 1; eps0 = 1;
    while(eps0 > 0.0000001 & xt<25){
      y1 = c0+xt*exp(xt)-c0*exp(xt)
      yt = (1+xt-c0)*exp(xt)
      xt_new = xt - y1/yt
      eps0 = abs(xt-xt_new)
      xt = xt_new
      k0 = k0+1
    }
    xt
  }
  sam_ztp_lam_fun = function(w_start,beta_start){
    lam = exp(w_start%*%beta_start  )
    lam = as.matrix(lam)
    x_new = apply(lam,1,sam_ztp)
  }
  sam_ztp = function(lam_k){
    sam_poi = rpois(1,lam_k)
    while(sam_poi==0){
      sam_poi = rpois(1,lam_k)
    }
    # sam_poi = sam_poi + rpois(1,0.3) #### add some noise
    sam_poi
  }
  sam_ztp_fun = function(w_start,beta_start){
    c_matrix = exp(w_start%*%beta_start )+1
    lam = apply(c_matrix,1,newton_lam)
    lam = as.matrix(lam)
    x_new = apply(lam,1,sam_ztp)
  }
  mean_ztp = function(lam_0){
    lam_0/(1-exp(-lam_0))
  }
  var_ztp = function(lam_0){
    lam_0/(1-exp(-lam_0)) - exp(-lam_0)*lam_0^2/(1-exp(-lam_0))^2
  }
  if(samp == 1){
    beta_ori = Beta_ini
  }else{
    beta_ori = R_mean0[ceiling(n/1000),]
  }
  
  X_ztp = sam_ztp_fun(W,beta_ori)
  
  dat = cbind(X_ztp,W[,-1])
  dat_orignal = as.data.frame(dat)
  n_row =nrow(dat_orignal);
  p =ncol(dat_orignal);
  beta = rep(0,p); k = 1; eps = 1;
  while(eps>0.00000001 && k<1000){
    c_mat = exp(W%*%beta)+1;
    lam = apply(c_mat,1,newton_lam);
    h1 = X_ztp/lam - 1/(1-exp(-lam));
    h2 = (1-exp(-lam)-lam*exp(-lam))/(1-exp(-lam))^2;
    h1_prime = exp(-lam)/(1-exp(-lam))^2 - X_ztp/lam^2
    h2_prime = exp(-lam)*(lam-2+(2+lam)*exp(-lam))/(1-exp(-lam))^3
    mu_new = c_mat%*%matrix(1,1,p)
    h3 = (mu_new-1)*W
    a1 = (h1_prime*h2-h2_prime*h1)/h2^3
    a2 = h1/h2
    G_matrix = matrix(0,p,p)
    partial_loglik = rep(0,p)
    for(ts in c(1:n_row)){
      G_matrix = G_matrix + a1[ts]*(h3[ts,]%*%t(h3[ts,])) + a2[ts]*(c_mat[ts]-1)*(W[ts,]%*%t(W[ts,]))
      partial_loglik = partial_loglik + a2[ts]*h3[ts,]
    }
    I_inf = ginv(G_matrix)
    beta = beta  -  I_inf%*%partial_loglik
    if(k>1){
      eps =  sum(abs(beta- beta_old))
    }
    beta_old = beta
    k = k+1
  }
  results0 =  c(beta,ceiling(n/1000))
}

p = 4
# Beta_zzz = matrix(c(0.2,0.2,0.5,0.5,0.5,-0.5),3,2,byrow = T)
Beta_zzz = matrix(c(0.2,0.2,0.5,0.5,0.2,-0.2,0.5,-0.5,-1,2,-1,2),3,p,byrow = T)

cv_pro = matrix(rep(0,3*3*p),3*3,p,byrow = T)
ciw = matrix(rep(0,3*3*p),3*3,p,byrow = T)

bis = matrix(rep(0,3*3*p),3*3,p,byrow = T)
mse = matrix(rep(0,3*3*p),3*3,p,byrow = T)

start_time = Sys.time()
numss = c(100,200,400)
for(zzz in c(1:3)){
  Beta_ini = Beta_zzz[zzz,]
  for(sss in c(1:3)){
    n1 = numss[sss]
    set.seed(100)
    W = cbind(matrix(1,n1,1),matrix(rnorm(n1,0.4,sqrt(0.1)),n1,1),
              matrix(floor(runif(n1,0,2)),n1,1),matrix((0.2+0.06*rt(n1,4)),n1,1))
    # W = cbind(matrix(1,n1,1),matrix(rnorm(n1,0.4,sqrt(0.1)),n1,1) )
    W = as.matrix(W)
    samp = 1
    ntimes = 10000
    no_cores <- detectCores()
    cl<-makeCluster(no_cores) 
    
    clusterEvalQ(cl, library(MASS))
    clusterExport(cl, "W");
    clusterExport(cl, "Beta_ini");
    clusterExport(cl, "samp");
    res_mean <- parLapply(cl, 1:ntimes, est_beta_fun)
    stopCluster(cl)
    R_mean0<- matrix(unlist(res_mean),ntimes,(p+1),byrow=T)
    R_mean0 = R_mean0[,c(1:p)]
    
    samp = 0
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    clusterEvalQ(cl, library(MASS))
    clusterExport(cl, "W");
    clusterExport(cl, "R_mean0");
    clusterExport(cl, "samp");
    res_mean1 <- parLapply(cl, 1:(ntimes*1000), est_beta_fun)
    stopCluster(cl)
    R_mean_all<- matrix(unlist(res_mean1),ntimes*1000,(p+1),byrow=T)
    
    bis[((zzz-1)*3+sss),] = apply(R_mean0,2,mean)-Beta_ini
    mse[((zzz-1)*3+sss),] = (apply(R_mean0,2,mean)-Beta_ini)^2 + apply(R_mean0,2,var)
    
    for(ttt in c(1:ntimes)){
      R_mean <- R_mean_all[R_mean_all[,(p+1)]==ttt,c(1:p)]
      ci = confident_interval_fun(R_mean)
      ciw[((zzz-1)*3+sss),] = ciw[((zzz-1)*3+sss),] +  ci[2,] - ci[1,]
      for(kkk in c(1:p)){
        if(ci[1,kkk]<=Beta_zzz[zzz,kkk] & ci[2,kkk]>=Beta_zzz[zzz,kkk]){
          cv_pro[((zzz-1)*3+sss),kkk] =  cv_pro[((zzz-1)*3+sss),kkk] + 1
        }
      }
    }
    print(c(zzz,sss))
  }
}
end_time = Sys.time()
print(end_time-start_time)


### accuracy of point estimates and interval estimates
mat_accuracy = matrix(0,16,9)
mat_accuracy[c(1:4)*4-3,] = t(bis)
mat_accuracy[c(1:4)*4-2,] = t(mse)
mat_accuracy[c(1:4)*4-1,] = t(cv_pro/ntimes)
mat_accuracy[c(1:4)*4,] = t(ciw/ntimes)


# rm(list = ls());gc()
