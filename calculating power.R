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
    while(eps0 > 0.00001 & xt<25){
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
  Inf.matrix = function(W,beta,X_ztp){
    X.ztp = X_ztp
    p = length(beta); n_row = dim(W)[1]
    c_mat = exp(W%*%beta)+1; mu.new = c_mat%*%matrix(1,1,p); lam.new = apply(c_mat,1,newton_lam);
    EX_ztp = apply(matrix(lam.new,n_row,1),1,mean_ztp);
    
    h1.new = X.ztp/lam.new - 1/(1-exp(-lam.new)); h1_prime.new = exp(-lam.new)/(1-exp(-lam.new))^2 - X.ztp/lam.new^2
    Eh1_prime.new = exp(-lam.new)/(1-exp(-lam.new))^2 - EX_ztp/lam.new^2;
    h2.new = (1-exp(-lam.new)-lam.new*exp(-lam.new))/(1-exp(-lam.new))^2; h2_prime.new = exp(-lam.new)*(lam.new-2+(2+lam.new)*exp(-lam.new))/(1-exp(-lam.new))^3
    h3.new = (mu.new-1)*W
    
    a1.new = (h1_prime.new*h2.new-h2_prime.new*h1.new)/h2.new^3; a2.new = h1.new/h2.new
    Ea1.new = Eh1_prime.new/h2.new^2
    
    H_matrix.new = matrix(0,p,p); I_matrix.new = matrix(0,p,p); partial_loglik = rep(0,p)
    for(ts in c(1:n_row)){
      I_matrix.new = I_matrix.new - Ea1.new[ts]*(h3.new[ts,]%*%t(h3.new[ts,]));
      H_matrix.new = H_matrix.new + a1.new[ts]*(h3.new[ts,]%*%t(h3.new[ts,])) + a2.new[ts]*(c_mat[ts]-1)*(W[ts,]%*%t(W[ts,]));
      partial_loglik = partial_loglik + a2.new[ts]*h3.new[ts,];
    }
    list(H_matrix.new,partial_loglik ,I_matrix.new)
  }
  Inf.matrix.lam = function(W,beta,X_ztp){
    H_matrix.new = matrix(0,p,p); I_matrix.new = matrix(0,p,p); partial_loglik = rep(0,p)
    lam.new = exp(W%*%beta);
    a1 = X_ztp/lam.new - 1 -1/(exp(lam.new)-1);
    a2 = -X_ztp/lam.new^2 + exp(lam.new)/(exp(lam.new)-1)^2;
    X_ztp.mean = apply(lam.new,1,mean_ztp)
    a3 = X_ztp.mean/lam.new - 1 -1/(exp(lam.new)-1);
    a4 = -X_ztp.mean/lam.new^2 + exp(lam.new)/(exp(lam.new)-1)^2;
    
    
    for(ts in c(1:n_row)){
      I_matrix.new = I_matrix.new - W[ts,]%*%t(W[ts,])*(a3[ts]*lam.new[ts]+a4[ts])*lam.new[ts];
      H_matrix.new = H_matrix.new + W[ts,]%*%t(W[ts,])*(a1[ts]*lam.new[ts]+a2[ts])*lam.new[ts];
      
      partial_loglik = partial_loglik + a1[ts]*lam.new[ts]*W[ts,];
    }
    list(H_matrix.new,partial_loglik ,I_matrix.new)
  }
  
  beta = Beta_ini
  
  if(mean.pwr == 1){ X_ztp = sam_ztp_fun(W,beta) }else{X_ztp = sam_ztp_lam_fun(W,beta)}
  dat = cbind(X_ztp,W[,-1])
  dat_orignal = as.data.frame(dat)
  n_row =nrow(dat_orignal)
  p =ncol(dat_orignal)
  tail.data = apply(as.matrix(X_ztp,nrow = n_row,byrow = F),1, function(x_z){sum(log(c(1:x_z)))})
  
  if(mean.pwr == 1){
    beta = rep(0,p); beta.old = rep(1,p); k = 0; eps = 1;# rr= 0.9;I_inf.old = matrix(0,p,p); G_matrix.old = matrix(0,p,p);
    
    while(eps>0.00000001 && k<1000){
      eps =  sum(abs(beta- beta.old))
      beta.old = beta; k = k+1;
      zzz = Inf.matrix(W,beta,X_ztp)
      H_matrix = matrix(unlist(zzz[1]),p,p,byrow=T)
      partial_loglik = matrix(unlist(zzz[2]),p,1,byrow=T)
      I_matrix = matrix(unlist(zzz[3]),p,p,byrow=T)
      inf.inv = ginv(H_matrix); 
      beta = beta  -  inf.inv%*%partial_loglik;
    }
    c_mat = exp(W%*%beta)+1;
    lam = apply(c_mat,1,newton_lam);
    
    lam.H01 = newton_lam(mean(X_ztp))
    beta.H01 = log(lam.H01/(1-exp(-lam.H01))-1)
    beta.H0 = c(beta.H01,rep(0,(p-1)))
    c_mat.H0 = exp(W%*%beta.H0)+1;
    lam.H0 = apply(c_mat.H0,1,newton_lam);
    
    ## hypothesis testing
    zzz.H0 = Inf.matrix(W,beta.H0,X_ztp)
    H_matrix.H0 = matrix(unlist(zzz.H0[1]),p,p,byrow=T)
    partial_loglik.H0 = matrix(unlist(zzz.H0[2]),p,1,byrow=T)
    I_matrix.H0 = matrix(unlist(zzz.H0[3]),p,p,byrow=T)
    
  }else{
    m1 <- vglm(X_ztp ~ ., family = pospoisson(), data = dat_orignal)
    beta = as.vector(t(coef(summary(m1))[, 1:2]))[c(1:p)*2-1]
    lam = exp(W%*%beta);
    
    lam.H01 = newton_lam(mean(X_ztp))
    beta.H01 = log(lam.H01)
    beta.H0 = c(beta.H01,rep(0,(p-1)))
    lam.H0 = exp(W%*%beta.H0)
    
    ## hypothesis testing
    zzz = Inf.matrix.lam(W,beta,X_ztp)
    H_matrix = matrix(unlist(zzz[1]),p,p,byrow=T)
    partial_loglik = matrix(unlist(zzz[2]),p,1,byrow=T)
    I_matrix = matrix(unlist(zzz[3]),p,p,byrow=T)
    
    zzz.H0 = Inf.matrix.lam(W,beta.H0,X_ztp)
    H_matrix.H0 = matrix(unlist(zzz.H0[1]),p,p,byrow=T)
    partial_loglik.H0 = matrix(unlist(zzz.H0[2]),p,1,byrow=T)
    I_matrix.H0 = matrix(unlist(zzz.H0[3]),p,p,byrow=T)
    
  }
  
  log_likelihood = sum(-log(1-exp(-lam))+X_ztp*log(lam) - lam - tail.data)
  log_likelihood.H0 = sum(-log(1-exp(-lam.H0))+X_ztp*log(lam.H0) - lam.H0 - tail.data)
  
  ## likelihood ratio test
  T.likelihood = 2*(log_likelihood - log_likelihood.H0)
  ## wald test
  T.Wald = t(beta-beta.H0)%*%I_matrix%*%(beta-beta.H0)
  ## score test
  T.score = t(partial_loglik.H0)%*%ginv(I_matrix.H0)%*%partial_loglik.H0
  
  # c_mat = exp(W%*%beta)+1;
  # lam_predict = apply(c_mat,1,newton_lam);
  # X_predict = sapply(lam_predict,mean_ztp,simplify = T)
  # mse = mean((X_predict-X_ztp)^2)
  
  aic = -2*log_likelihood + 2*p
  bic = p*log(n_row) - 2*log_likelihood
  
  results0 =  c(beta,T.likelihood, T.Wald,T.score,aic,bic,ceiling(n/10000))
  # names(results0) = c(rep("beta",p),"likelihoodratio","wald","score","aic","bic")
  results0
}

library(parallel)
# Beta_zzz = matrix(c(0,0,0.2,0.2,0.4,0.4),3,2,byrow = T)
# cv_pro = matrix(rep(0,3*3*8),3*8,3,byrow = T)
# p = 2
cv_pro = c()
Beta_zzz = matrix(c(0,0,0,0,
                    -0.5,-0.2, 0.4, -0.2,
                    0.3,0.2, 0.2, -0.2),3,4,byrow = T)

p = 4
start_time = Sys.time()
numss = c(50,100,150,200,250,300,350,400)
stands = qchisq(c(0.025,0.975),p-1)  ## chi-square value
for(zzz.count in c(1:3)){
  Beta_ini = Beta_zzz[zzz.count,]
  for(sss.count in c(1:8)){
    n1 = numss[sss.count]
    set.seed(100)
    # W = cbind(matrix(1,n1,1),matrix(runif(n1),n1,1))
    
    W = cbind(matrix(1,n1,1),matrix(rnorm(n1,0.5,sqrt(0.1)),n1,1))
    W = cbind(W, matrix(floor(runif(n1,0,2)),n1,1))
    W = cbind(W, 0.5+0.6*matrix(rt(n1,5),n1,1))
    W = as.matrix(W)
    
    ntimes = 1
    mtimes = 10000
    mean.pwr = 0
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    clusterEvalQ(cl, library(MASS))
    clusterEvalQ(cl, library(VGAM));
    clusterExport(cl, "mean.pwr");
    clusterExport(cl, "W");
    clusterExport(cl, "Beta_ini");
    clusterExport(cl, "ntimes");
    res_mean1 <- parLapply(cl, 1:(ntimes*mtimes), est_beta_fun)
    stopCluster(cl)
    R_mean_all<- matrix(unlist(res_mean1),ntimes*mtimes,(p+6),byrow=T)
    
    R_mean <- R_mean_all[,c((p+1):(p+3))]
    r_k = colSums(R_mean>stands[2])  + colSums(R_mean<stands[1])
    if(zzz.count==1 & sss.count ==1){cv_pro = c(Beta_ini,n1,r_k/mtimes) } else{cv_pro = cbind(cv_pro,c(Beta_ini,n1,r_k/mtimes))}
    
  }
}
end_time = Sys.time()
print(end_time-start_time)
cv_pro
cv_pro_lam = cv_pro


for(zzz.count in c(1:3)){ 
  
  Beta_ini = Beta_zzz[zzz.count,]
  for(sss.count in c(1:8)){
    n1 = numss[sss.count]
    set.seed(100)
    # W = cbind(matrix(1,n1,1),matrix(runif(n1),n1,1))
    
    W = cbind(matrix(1,n1,1),matrix(rnorm(n1,0.5,sqrt(0.1)),n1,1))
    W = cbind(W, matrix(floor(runif(n1,0,2)),n1,1))
    W = cbind(W, 0.5+0.6*matrix(rt(n1,5),n1,1))
    W = as.matrix(W)
    
    ntimes = 1
    mtimes = 10000
    mean.pwr = 1
    no_cores <- detectCores()
    cl<-makeCluster(no_cores)
    clusterEvalQ(cl, library(MASS))
    clusterEvalQ(cl, library(VGAM));
    clusterExport(cl, "mean.pwr");
    clusterExport(cl, "W");
    clusterExport(cl, "Beta_ini");
    clusterExport(cl, "ntimes");
    res_mean1 <- parLapply(cl, 1:(ntimes*mtimes), est_beta_fun)
    stopCluster(cl)
    R_mean_all<- matrix(unlist(res_mean1),ntimes*mtimes,(p+6),byrow=T)
    
    R_mean <- R_mean_all[,c((p+1):(p+3))]
    r_k = colSums(R_mean>stands[2])  + colSums(R_mean<stands[1])
    if(zzz.count==1 & sss.count ==1){cv_pro = c(Beta_ini,n1,r_k/mtimes) } else{cv_pro = cbind(cv_pro,c(Beta_ini,n1,r_k/mtimes))}
    
  }
}
end_time = Sys.time()
print(end_time-start_time)
cv_pro
cv_pro_mu = cv_pro

# ### plot power
# cv_pro = cv_pro_mu
plot(numss[c(1:8)],cv_pro[6, c(1:8)], ylim = c(0,0.1),type = "l",xlab = "data size", ylab="empirical significant level")
lines(numss[c(1:8)],cv_pro[7, c(1:8)],col="blue",lty=2)
lines(numss[c(1:8)],cv_pro[8, c(1:8)],col="red",lty=3)
legend('bottomright',inset = .1,expression("likelihood ratio","Wald","score"),lty=c(1,2,3),col = c("black","blue","red"))


plot(numss[c(1:8)],cv_pro[6, c(1:8)+8], ylim = c(0,1),type = "l",xlab = "data size", ylab="empirical power")
lines(numss[c(1:8)],cv_pro[7, c(1:8)+8],col="blue",lty=2)
lines(numss[c(1:8)],cv_pro[8, c(1:8)+8],col="red",lty=3)
legend('bottomright',inset = .1,expression("likelihood ratio","Wald","score"),lty=c(1,2,3),col = c("black","blue","red"))

plot(numss[c(1:8)],cv_pro[6, c(1:8)+16], ylim = c(0,1),type = "l",xlab = "data size", ylab="empirical power")
lines(numss[c(1:8)],cv_pro[7, c(1:8)+16],col="blue",lty=2)
lines(numss[c(1:8)],cv_pro[8, c(1:8)+16],col="red",lty=3)
legend('bottomright',inset = .1,expression("likelihood ratio","Wald","score"),lty=c(1,2,3),col = c("black","blue","red"))
