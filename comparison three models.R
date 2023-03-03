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
  sam_ztp_lam_fun = function(w_start,beta_start){
    lam = exp(w_start%*%beta_start  )
    lam = as.matrix(lam)
    x_new = apply(lam,1,sam_ztp)
  }
  sam_ztp_mu_fun = function(w_start,beta_start){
    c_matrix = exp(w_start%*%beta_start)
    lam = apply(c_matrix,1,newton_lam_k)
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
    lam = apply(c_matrix,1,newton_lam_k)
    lam = as.matrix(lam)
    x_new = apply(lam,1,sam_ztp)
  }
  
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
    beta = c(log(trun.k)+1,rep(0,(p-1))); step_len = 0.001; k = 1; eps = 1; log_likelihood = 0;
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
        # beta.new = beta  - step.len*I_inf%*%partial_loglik
        beta.new = beta  + step.len*partial_loglik
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
      beta.new = beta - I_inf%*%partial_loglik
      c_mat = exp(W%*%beta.new);
      lam = c_mat;
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

  
  trun.k = 1 ## zero truncated poisson
  
  simu = 1 ##  simu = 1, if simulation; otherwise true data;
  interception = 1 ##  interception = 1, if consider interception;
  mean.pwr = 0 ##  mean.pwr = 1, calculate log mean based regression; otherwise, log lambda based regression;
  log_mu.pwr = 1
  log_lambda.pwr = 0
  markbreak = 0
  
  dat_orignal = W  # original data
  n_row =nrow(dat_orignal)
  
  if(interception == 1){
    p =ncol(dat_orignal); W = cbind(matrix(1,n_row,1),dat_orignal[,-1]); W = as.matrix(W)
  }else{
    p =ncol(dat_orignal)-1; W = dat_orignal[,-1]; W = as.matrix(W)
  }
  
  
  if(simu==1){
    beta = Beta_ini
    if(mean.pwr == 1){
      X_ztp = sam_ztp_fun(W,beta)
    }else if(log_mu.pwr == 1){
      X_ztp = sam_ztp_mu_fun(W,beta)
    }else if(log_lambda.pwr == 1){
      X_ztp = sam_ztp_lam_fun(W,beta)
    }
    # dat = cbind(X_ztp,W[,-1])
    # for(i in c(1:p)){
    #   mark1 = sum(dat[dat[,i]>0,1])/sum(dat[,i]>0)
    #   if(mark1==1){dat[dat[,i]>0,1]=dat[dat[,i]>0,1]+0.01}
    # }
    # dat_orignal = as.data.frame(dat); X_ztp = dat_orignal[,1];
  }else{
    X_ztp = dat_orignal[,1];
    dat = cbind(X_ztp,W[,-1])
    dat_orignal = as.data.frame(dat)
  }
  
  
  tail.data = apply(as.matrix(X_ztp,nrow = n_row,byrow = F),1, function(x_z){sum(log(c(1:x_z)))})
  C = diag(1,100,100)[-1,] - rbind(diag(1,100,100)[-c(1,2),],rep(0,100))
  C = C[c(1:(p-1)),c(1:p)]
  
  
  results.log.mu = est_ztp_log_mu_fun(W,X_ztp,trun.k)
  results.log.mu_k = est_ztp_log_mu_k_fun(W,X_ztp,trun.k)
  results.log.lambda = est_ztp_log_lambda_fun(W,X_ztp,trun.k)
  
  results.case = c(results.log.mu,
                       results.log.mu_k,
                       results.log.lambda)
  
  
}

# rm(list=ls());gc()

set.seed(100)
n1 = 1000
p = 3
# W = cbind(matrix(1,n1,1),matrix(rnorm(n1*(p-1),0.4,sqrt(0.5)),n1,(p-2)),matrix(floor(runif(n1,0,2)),n1,1))
W = cbind(matrix(1,n1,1),matrix(runif(n1*(p-1)),n1,(p-2)),matrix(floor(runif(n1,0,2)),n1,1))
W = as.matrix(W)

Beta_ini = c(0.5,1,0.3)

library(parallel)
no_cores <- detectCores()
cl<-makeCluster(no_cores)
clusterEvalQ(cl, library(MASS))
clusterEvalQ(cl, library(VGAM))
clusterExport(cl, "W");
clusterExport(cl, "Beta_ini");
ntimes = 10000
start_time = Sys.time()
res_mean_sigle <- parLapply(cl, 1:ntimes,  est_beta_fun)
end_time = Sys.time()
print(end_time-start_time)
stopCluster(cl)


# res_mean_sigle.old = res_mean_sigle

# estimation of parameters for two models
R_mean_sigle<- matrix(unlist(res_mean_sigle),ntimes,(p*4+7)*3,byrow=T)
ci = confident_interval_fun(R_mean_sigle)
mea = apply(R_mean_sigle,2,mean)
sd1 = apply(R_mean_sigle,2,sd)

comparison2 = rbind(ci,mea,sd1)
write.csv(comparison2,'C:\\Users\\Lenovo\\Desktop\\comparison3.csv')

