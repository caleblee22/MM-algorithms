##########################################################
####### snake data set ###################################
##########################################################
dat <- read.table("C:\\Users\\Caleb Lee\\Desktop\\papers\\ZTP_reg final\\ztp regression\\Snakes.txt",header = T)
# dat1 = dat[,c(5,7,9,10)]
library("nnet")
dat.season = class.ind(dat$Season)
dat[,"Road_Loc"] = c(rep(0,83),rep(1,47))
dat = cbind.data.frame(dat,dat.season)
trun.k = 1
dat1 = dat[c(1:130)[dat$N_days>=trun.k],c(5,7,10,13,14)]
# dat1$PDayRainTot_Rain <- dat1$PDayRain*dat1$Tot_Rain
p = dim(dat1)[2]
W = dat1[,c(2:p)]
w1 = matrix(1,nrow = dim(W)[1], 1)
W = cbind(w1,W)
W = as.matrix(W)
dat_orignal = as.data.frame(dat1)
n_row =nrow(dat_orignal);
p =ncol(dat_orignal);
X_ztp = dat_orignal[,1]
names(dat_orignal)[1] = 'X_ztp'


##########################################################
####### snake data set ###################################
##########################################################
dat <- read.csv("C:\\Users\\Caleb Lee\\Desktop\\papers\\ZTP_reg final\\ztp regression\\SEER Breast Cancer Dataset.csv",header = T)
dat = dat[dat$Status=='Dead',]
dat = dat[,c(15,1,3,9,11,12)]
library("nnet")
Marital = class.ind(dat$Marital.Status)[,2]
A.Stage = class.ind(dat$A.Stage)[,2]
Progesterone = class.ind(dat$Progesterone.Status)[,2]
Estrogen = class.ind(dat$Estrogen.Status)[,2]

dat1 = cbind( ceiling(dat[,1]/12),dat[,2],A.Stage,Estrogen,Progesterone,Marital )

p = dim(dat1)[2]
W = dat1[,c(2:p)]
w1 = matrix(1,nrow = dim(W)[1], 1)
W = cbind(w1,W)
W = as.matrix(W)
dat_orignal = as.data.frame(dat1)
n_row =nrow(dat_orignal);
p =ncol(dat_orignal);
X_ztp = dat_orignal[,1]
names(dat_orignal)[1] = 'X_ztp'



factor.xxx = function(xxx){
  sum(log(c(1:xxx)))
}
newton_lam_k = function(c0){
  if(c0<0){c0=c0+1}
  xt = c0; k0 = 1; eps0 = 1;
  if(trun.k==1){
    while(eps0 > 0.0000001 & xt<25){
      y1 = c0+xt*exp(xt)-c0*exp(xt)
      yt = (1+xt-c0)*exp(xt)
      xt_new = xt - y1/yt
      eps0 = abs(xt-xt_new)
      xt = xt_new
      k0 = k0+1
    }
  }else{
    while(eps0 > 0.0000001 & xt<25){
      tra.data = xt^c(0:(trun.k-1))*exp(-xt)/factorial(c(0:(trun.k-1)))
      xt_new = c0*(1-sum(tra.data))/(1-sum(tra.data) + tra.data[trun.k])
      eps0 = abs(xt-xt_new)
      xt = xt_new
      k0 = k0+1
      xt
    }
  }
  xt
}
log_likelihood_fun = function(X_ztp,lam,trun.k){
  n_row = length(X_ztp)
  tra.data = apply(matrix(lam,length(lam),1), 1, function(lam1) if(lam1<40){ sum(lam1^c(0:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-1))))}else{0})
  log_likelihood_old = sum(-log(1-tra.data)+X_ztp*log(lam) - lam - apply(matrix(X_ztp,n_row,1),1,factor.xxx))
}
Inf.matrix_log_mu_trun.k = function(W,beta,X_ztp,trun.k){
  X.ztp = X_ztp
  p = length(beta); n_row = dim(W)[1]
  c_mat = exp(W%*%beta); mu.new = c_mat%*%matrix(1,1,p); lam.new = apply(c_mat,1,newton_lam_k);
  A1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-1)))))
  A2 = exp(-lam.new)*lam.new^(trun.k-1)/factorial(trun.k-1)
  if((trun.k-1)>0){
    A3 = exp(-lam.new)*(lam.new^(trun.k-2)/factorial(trun.k-2) - lam.new^(trun.k-1)/factorial(trun.k-1))
  }else{
    A3 = -exp(-lam.new)
  }
  if((trun.k-1)>0){
    B1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) lam1 - sum(lam1^c(1:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-2)))))
  }else{
    B1 = lam.new
  }
  if((trun.k-1)>0){
    B2 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-2))*exp(-lam1)/factorial(c(0:(trun.k-2)))) +
                 lam1^(trun.k-1)*exp(-lam1)/factorial(trun.k-2))
  }else{
    B2 = 1
  }
  if((trun.k-1)>0){
    B3 = apply(matrix(lam.new,n_row,1), 1, function(lam1)  (-lam1+trun.k)*exp(-lam1)*lam1^c(trun.k-2)/factorial(trun.k-2))
  }else{
    B3 =0
  }
  h1.new = X.ztp/lam.new - 1-A2/A1
  h1_prime.new = -(A3*A1-A2^2)/A1^2 - X.ztp/lam.new^2
  h2.new = (B2*A1-B1*A2)/A1^2
  h2_prime.new = ((B3*A1-B1*A3)*A1 - 2*A2*(B2*A1-B1*A2))/A1^3
  h3.new = mu.new*W
  a1.new = (h1_prime.new*h2.new-h2_prime.new*h1.new)/h2.new^3; a2.new = h1.new/h2.new
  H_matrix.new = matrix(0,p,p); partial_loglik = rep(0,p)
  # h2.new.old = h2.new
  # lam.new = lam.new+0.000001
  # (h2.new - h2.new.old)/0.000001
  for(ts in c(1:n_row)){
    H_matrix.new = H_matrix.new + a1.new[ts]*(h3.new[ts,]%*%t(h3.new[ts,])) + a2.new[ts]*c_mat[ts]*(W[ts,]%*%t(W[ts,]));
    partial_loglik = partial_loglik + a2.new[ts]*h3.new[ts,];
  }
  list(H_matrix.new,partial_loglik)
}
Inf.matrix_trun.k = function(W,beta,X_ztp,trun.k){
  X.ztp = X_ztp
  p = length(beta); n_row = dim(W)[1]
  c_mat = exp(W%*%beta)+trun.k; mu.new = c_mat%*%matrix(1,1,p); lam.new = apply(c_mat,1,newton_lam_k);
  
  A1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-1)))))
  A2 = exp(-lam.new)*lam.new^(trun.k-1)/factorial(trun.k-1)
  if((trun.k-1)>0){
    A3 = exp(-lam.new)*(lam.new^(trun.k-2)/factorial(trun.k-2) - lam.new^(trun.k-1)/factorial(trun.k-1))
  }else{
    A3 = -exp(-lam.new)
  }
  if((trun.k-1)>0){
    B1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) lam1 - sum(lam1^c(1:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-2)))))
  }else{
    B1 = lam.new
  }
  if((trun.k-1)>0){
    B2 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-2))*exp(-lam1)/factorial(c(0:(trun.k-2)))) +
                 lam1^(trun.k-1)*exp(-lam1)/factorial(trun.k-2))
  }else{
    B2 = 1
  }
  if((trun.k-1)>0){
    B3 = apply(matrix(lam.new,n_row,1), 1, function(lam1)  (-lam1+trun.k)*exp(-lam1)*lam1^c(trun.k-2)/factorial(trun.k-2))
  }else{
    B3 =0
  }
  h1.new = X.ztp/lam.new - 1-A2/A1
  h1_prime.new = -(A3*A1-A2^2)/A1^2 - X.ztp/lam.new^2
  h2.new = (B2*A1-B1*A2)/A1^2
  h2_prime.new = ((B3*A1-B1*A3)*A1 - 2*A2*(B2*A1-B1*A2))/A1^3
  h3.new = (mu.new-trun.k)*W
  a1.new = (h1_prime.new*h2.new-h2_prime.new*h1.new)/h2.new^3; a2.new = h1.new/h2.new
  H_matrix.new = matrix(0,p,p); partial_loglik = rep(0,p)
  # h2.new.old = h2.new
  # lam.new = lam.new+0.000001
  # (h2.new - h2.new.old)/0.000001
  for(ts in c(1:n_row)){
    H_matrix.new = H_matrix.new + a1.new[ts]*(h3.new[ts,]%*%t(h3.new[ts,])) + a2.new[ts]*(c_mat[ts]-trun.k)*(W[ts,]%*%t(W[ts,]));
    partial_loglik = partial_loglik + a2.new[ts]*h3.new[ts,];
  }
  list(H_matrix.new,partial_loglik)
}
Inf.matrix_log_lambda_trun.k = function(W,beta,X_ztp,trun.k){
  X.ztp = X_ztp
  p = length(beta); n_row = dim(W)[1]
  c_mat = exp(W%*%beta);
  mu.new = c_mat%*%matrix(1,1,p);
  lam.new = c_mat;
  
  
  A1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-1)))))
  A2 = exp(-lam.new)*lam.new^(trun.k-1)/factorial(trun.k-1)
  
  if((trun.k-1)>0){
    A3 = exp(-lam.new)*(lam.new^(trun.k-2)/factorial(trun.k-2) - lam.new^(trun.k-1)/factorial(trun.k-1))
  }else{
    A3 = -exp(-lam.new)
  }
  if((trun.k-1)>0){
    B1 = apply(matrix(lam.new,n_row,1), 1, function(lam1) lam1 - sum(lam1^c(1:(trun.k-1))*exp(-lam1)/factorial(c(0:(trun.k-2)))))
  }else{
    B1 = lam.new
  }
  if((trun.k-1)>0){
    B2 = apply(matrix(lam.new,n_row,1), 1, function(lam1) 1 - sum(lam1^c(0:(trun.k-2))*exp(-lam1)/factorial(c(0:(trun.k-2)))) +
                 lam1^(trun.k-1)*exp(-lam1)/factorial(trun.k-2))
  }else{
    B2 = 1
  }
  if((trun.k-1)>0){
    B3 = apply(matrix(lam.new,n_row,1), 1, function(lam1)  (-lam1+trun.k)*exp(-lam1)*lam1^c(trun.k-2)/factorial(trun.k-2))
  }else{
    B3 =0
  }
  h1.new = X.ztp/lam.new - 1-A2/A1
  h1_prime.new = -(A3*A1-A2^2)/A1^2 - X.ztp/lam.new^2
  h3.new = mu.new*W
  a1.new = h1_prime.new
  a2.new = h1.new
  H_matrix.new = matrix(0,p,p); partial_loglik = rep(0,p)
  for(ts in c(1:n_row)){
    H_matrix.new = H_matrix.new + a1.new[ts]*(h3.new[ts,]%*%t(h3.new[ts,])) + a2.new[ts]*c_mat[ts]*(W[ts,]%*%t(W[ts,]));
    partial_loglik = partial_loglik + a2.new[ts]*h3.new[ts,];
  }
  list(H_matrix.new,partial_loglik)
}
mean_ztp_trun.k = function(lam_0,trun.k=1){
  if(lam_0<40){
    tra.data = lam_0^c(0:(trun.k-1))*exp(-lam_0)/factorial(c(0:(trun.k-1)))
    if(trun.k>1){
      tra.data1 = tra.data[c(1:(trun.k-1))]
    }else{
      tra.data1 = 0
    }
    if(trun.k>2){
      tra.data2 = tra.data[c(1:(trun.k-2))]
    }else{
      tra.data2 = 0
    }
  }else{
    tra.data = 0
    tra.data1 = 0
    tra.data2 = 0
  }
  (lam_0 - sum(tra.data1))/(1-sum(tra.data))
}
var_ztp_trun.k = function(lam_0,trun.k=1){
  if(lam_0<40){
    tra.data = lam_0^c(0:(trun.k-1))*exp(-lam_0)/factorial(c(0:(trun.k-1)))
    if(trun.k>1){
      tra.data1 = tra.data[c(1:(trun.k-1))]
    }else{
      tra.data1 = 0
    }
    if(trun.k>2){
      tra.data2 = tra.data[c(1:(trun.k-2))]
    }else{
      tra.data2 = 0
    }
  }else{
    tra.data = 0
    tra.data1 = 0
    tra.data2 = 0
  }
  EX = (lam_0 - sum(tra.data1))/(1-sum(tra.data))
  (lam_0^2 - sum(tra.data2))/(1-sum(tra.data)) + EX*(1-EX)
}

est_ztp_log_mu_fun = function(W,X_ztp,trun.k){
  p = dim(W)[2]
  start_time = Sys.time()
  beta = c(log(trun.k)+1,rep(0,(p-1))); k = 1; eps = 1; log_likelihood = 0;
  step.len = 1
  c_mat = exp(W%*%beta);
  lam = apply(c_mat,1,newton_lam_k);
  log_likelihood_old = log_likelihood_fun(X_ztp,lam,trun.k)
  while(eps>1e-12){
    zzz = Inf.matrix_log_mu_trun.k(W,beta,X_ztp,trun.k)
    G_matrix = matrix(unlist(zzz[1]),p,p,byrow=T)
    partial_loglik = matrix(unlist(zzz[2]),p,1,byrow=T)
    I_inf = ginv(G_matrix)
    yes.no = TRUE
    step.len = 0.1
    while(yes.no){
      beta.new = beta  - step.len*I_inf%*%partial_loglik
      # beta.new = beta  + step.len*partial_loglik
      if((sum(W%*%beta.new < log(trun.k)) + sum(W%*%beta.new>100))>0){
        step.len = 0.5*step.len
      }else{
        c_mat = exp(W%*%beta.new);
        lam = apply(c_mat,1,newton_lam_k);
        log_likelihood = log_likelihood_fun(X_ztp,lam,trun.k)
        if(log_likelihood<log_likelihood_old){
          step.len = 0.5*step.len
        }else{
          yes.no = FALSE
        }
      }
    }
    if(k>1){
      eps =  abs(log_likelihood_old-log_likelihood)/(abs(log_likelihood_old)+1)
    }
    log_likelihood_old = log_likelihood
    beta = beta.new
    k = k+1
  }
  end_time = Sys.time()
  times = difftime(end_time,start_time,units='secs')
  
  X_predict = sapply(lam,mean_ztp_trun.k,simplify = T)
  X_predict_var = sapply(lam,var_ztp_trun.k,simplify = T)
  goodness_of_fit = sum((X_predict-X_ztp)^2/X_predict_var)
  aic = -2*log_likelihood + 2*p
  bic = -2*log_likelihood + p*log(length(lam))
  Wald.stat = t(beta)%*%(-G_matrix)%*%beta
  T.Wald.parameters = rep(0,p)
  p.parameters = rep(0,p)
  se.parameters = rep(0,p)
  for(i in c(1:p)){
    C = matrix(0,1,p)
    C[i] = 1
    T.Wald.parameters[i] = sqrt(t(C%*%beta)%*%ginv(C%*%ginv((-G_matrix))%*%t(C))%*%(C%*%beta))
    se.parameters[i] = sqrt(C%*%ginv((-G_matrix))%*%t(C))
    p.parameters[i] = 2*(1-pnorm(abs(beta[i]),0,se.parameters[i]))
  }
  
  
  results.log.mu =  c(beta,se.parameters,T.Wald.parameters,p.parameters,log_likelihood,aic,bic,Wald.stat,goodness_of_fit,k,times)
  # p.parameters
}
est_ztp_log_mu_k_fun = function(W,X_ztp,trun.k){
  p = dim(W)[2]
  start_time = Sys.time()
  beta = rep(0,p); step_len = 1; k = 1; eps = 1; log_likelihood = 0;
  c_mat = exp(W%*%beta)+trun.k;
  lam = apply(c_mat,1,newton_lam_k);
  log_likelihood_old = log_likelihood_fun(X_ztp,lam,trun.k)
  while(eps>1e-12){
    zzz = Inf.matrix_trun.k(W,beta,X_ztp,trun.k)
    G_matrix = matrix(unlist(zzz[1]),p,p,byrow=T)
    partial_loglik = matrix(unlist(zzz[2]),p,1,byrow=T)
    I_inf = ginv(G_matrix)
    beta.new = beta - I_inf%*%partial_loglik
    c_mat = exp(W%*%beta.new)+trun.k;
    lam = apply(c_mat,1,newton_lam_k);
    log_likelihood = log_likelihood_fun(X_ztp,lam,trun.k)
    eps =  abs(log_likelihood_old-log_likelihood)/(abs(log_likelihood_old)+1)
    log_likelihood_old = log_likelihood
    beta = beta.new
    k = k+1
  }
  end_time = Sys.time()
  times = difftime(end_time,start_time,units='secs')
  
  X_predict = sapply(lam,mean_ztp_trun.k,simplify = T)
  X_predict_var = sapply(lam,var_ztp_trun.k,simplify = T)
  goodness_of_fit = sum((X_predict-X_ztp)^2/X_predict_var)
  aic = -2*log_likelihood + 2*p
  bic = -2*log_likelihood + p*log(length(lam))
  Wald.stat = t(beta)%*%(-G_matrix)%*%beta
  T.Wald.parameters = rep(0,p)
  p.parameters = rep(0,p)
  se.parameters = rep(0,p)
  for(i in c(1:p)){
    C = matrix(0,1,p)
    C[i] = 1
    T.Wald.parameters[i] = sqrt(t(C%*%beta)%*%ginv(C%*%ginv((-G_matrix))%*%t(C))%*%(C%*%beta))
    se.parameters[i] = sqrt(C%*%ginv((-G_matrix))%*%t(C))
    p.parameters[i] = 2*(1-pnorm(abs(beta[i]),0,se.parameters[i]))
  }
  
  results.log.mu_k =  c(beta,se.parameters,T.Wald.parameters,p.parameters,log_likelihood,aic,bic,Wald.stat,goodness_of_fit,k,times)
  # p.parameters
}
est_ztp_log_lambda_fun = function(W,X_ztp,trun.k){
  p = dim(W)[2]
  start_time = Sys.time()
  
  beta = rep(0,p); step_len = 1; k = 1; eps = 1; log_likelihood = 0;
  c_mat = exp(W%*%beta);
  lam = c_mat;
  log_likelihood_old = log_likelihood_fun(X_ztp,lam,trun.k)
  while(eps>1e-12){
    zzz = Inf.matrix_log_lambda_trun.k(W,beta,X_ztp,trun.k)
    G_matrix = matrix(unlist(zzz[1]),p,p,byrow=T)
    partial_loglik = matrix(unlist(zzz[2]),p,1,byrow=T)
    I_inf = ginv(G_matrix)
    yes.no = TRUE
    step.len = 1
    while(yes.no){
      beta.new = beta  - step.len*I_inf%*%partial_loglik
      # beta.new = beta  + step.len*partial_loglik
      if((sum(W%*%beta.new< -5)+sum(W%*%beta.new>100))>0){
        step.len = 0.5*step.len
      }else{
        lam = exp(W%*%beta.new);
        log_likelihood = log_likelihood_fun(X_ztp,lam,trun.k)
        if(log_likelihood<log_likelihood_old){
          step.len = 0.5*step.len
        }else{
          yes.no = FALSE
        }
      }
    }
    # beta.new = beta - I_inf%*%partial_loglik
    lam = exp(W%*%beta.new);
    log_likelihood = log_likelihood_fun(X_ztp,lam,trun.k)
    eps =  abs(log_likelihood_old-log_likelihood)/(abs(log_likelihood_old)+1)
    log_likelihood_old = log_likelihood
    beta = beta.new
    k = k+1
  }
  end_time = Sys.time()
  times = difftime(end_time,start_time,units='secs')
  
  X_predict = sapply(lam,mean_ztp_trun.k,simplify = T)
  X_predict_var = sapply(lam,var_ztp_trun.k,simplify = T)
  goodness_of_fit = sum((X_predict-X_ztp)^2/X_predict_var)
  aic = -2*log_likelihood + 2*p
  bic = -2*log_likelihood + p*log(length(lam))
  Wald.stat = t(beta)%*%(-G_matrix)%*%beta
  T.Wald.parameters = rep(0,p)
  p.parameters = rep(0,p)
  se.parameters = rep(0,p)
  for(i in c(1:p)){
    C = matrix(0,1,p)
    C[i] = 1
    T.Wald.parameters[i] = sqrt(t(C%*%beta)%*%ginv(C%*%ginv((-G_matrix))%*%t(C))%*%(C%*%beta))
    se.parameters[i] = sqrt(C%*%ginv((-G_matrix))%*%t(C))
    p.parameters[i] = 2*(1-pnorm(abs(beta[i]),0,se.parameters[i]))
  }
  results.log.lambda =  c(beta,se.parameters,T.Wald.parameters,p.parameters,log_likelihood,aic,bic,Wald.stat,goodness_of_fit,k,times)
  # p.parameters
}


# names(dat1.revise.pre)
# W = dat1.revise.pre[dat1$Status=="Dead",-c(2,8,9,10,11,12,13,14,15,16)]
# W = as.matrix(W)
results.log.mu = est_ztp_log_mu_fun(W,X_ztp,trun.k)
results.log.mu_k = est_ztp_log_mu_k_fun(W,X_ztp,trun.k)
results.log.lambda = est_ztp_log_lambda_fun(W,X_ztp,trun.k)


# m1 <- vglm(X_ztp ~ ., family = pospoisson(), data = dat_orignal)

# write.csv(results.case2,"C:\\Users\\Caleb Lee\\Desktop\\papers\\Í¶¸å\\CSDA\\major revison\\results.case2.csv")
results.case = rbind(results.log.mu,
      results.log.mu_k,
      results.log.lambda)
results.case



