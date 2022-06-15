data {
    int n;//total number of subjects 
    int n_RE;// total numer of random effects
    int N1;// total number of observations 
    int ncx1;// total number of variables in the fixed part
    int id1[N1];// vector of repeated ids
    int RE_ind1[2];// vector with 2 elements
    vector[N1] y1;// the dependent variable 
    matrix[N1, 2] Z_1;// the design matrix of the random effects when we have raw variables 
    matrix[n, ncx1] Xhc1;// the design matrix of the fixed part when we have raw variables 
    real<lower=0> scale_betas1;// real number. This will be the standard error of the prior distribution of the coefficients 
    real<lower=0> scale_sigmas;// real number. This will be the standard error of the prior distribution of the standard deviation of the error terms  
    real<lower=0> scale_diag_D; // real number. This will be the standard error of the prior distribution of the standard deviation of the random terms 
    real<lower=0> lkj_shape;//real number. This will be the scale of the cholesky distribution
}
 
parameters {
    vector[ncx1] betas1;// coefficients of the variables in the fixed part 
    real<lower = 0> sigma1;// standard deviation of the error terms 
    matrix[n, n_RE] u;// random terms u 
    vector<lower = 0>[n_RE] L_var_D;// vector with n_RE elements. It will contain the variances of the random effects u
    cholesky_factor_corr[n_RE] L_corr_D;// cholesky factor correlation matrix
}
 
transformed parameters {
    vector[N1] eta1;// it will include the expected values of the model 
    matrix[n, n_RE] mu_u;
    for (i in 1:n) {
        mu_u[i, 1] = Xhc1[i, 1] * betas1[1] + Xhc1[i, 3] * betas1[3] + Xhc1[i, 4] * betas1[4];// it will give the mean
        // of the multivariate normal distribution of the random intercept of the u 
        mu_u[i, 2] = Xhc1[i, 2] * betas1[2];// it will give the mean
        // of the multivariate normal distribution of the random slope of the u
    }
    for (j1 in 1:N1) {
        eta1[j1] = Z_1[j1, 1] * u[id1[j1], 1] + Z_1[j1, 2] * u[id1[j1], 2];// estimated values 
    }
}

  // lets define now the prior distributions 
  
model {
    matrix[n_RE, n_RE] L_D;
    L_D = diag_pre_multiply(L_var_D, L_corr_D);
    L_var_D ~ cauchy(0, scale_diag_D);
    L_corr_D ~ lkj_corr_cholesky(lkj_shape);
    for (i in 1:n) {
        u[i, ] ~ multi_normal_cholesky(mu_u[i, ], L_D);
    }
    for (k1 in 1:ncx1) {
        betas1[k1] ~ normal(0.0, scale_betas1);
    }
     sigma1 ~ cauchy(0, scale_sigmas);
    y1 ~ normal(eta1, sigma1);
}
 
generated quantities {
    matrix[n_RE, n_RE] D;// variance covariance matrix of the random effects b
    matrix[n, n_RE] b;// these will be the random effects b 
    D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';// cholesky decomposition. We calculate the 
    //variance-covariance matrix of the random terms u 
    b = u - mu_u;// we calcuate the random terms b
}
