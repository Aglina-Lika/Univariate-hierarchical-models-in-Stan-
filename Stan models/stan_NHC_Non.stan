data {
    int n;//total number of subjects 
    int n_RE;// total numer of random effects
    int N1;// total number of observations
    int ncx1;// total number of variables in the fixed part
    int id1[N1];// vector of repeated ids
    int RE_ind1[2];// vector with 2 elements
    vector[N1] y1;// the dependent variable 
    matrix[N1, 2] Z_1;// the design matrix of the random part when we have raw variables 
    matrix[N1, ncx1] X1;// the design matrix of the fixed part when we have raw variables 
    real<lower=0> scale_betas1;// real number. This will be the standard error of the prior distribution of the coefficients 
    real<lower=0> scale_sigmas;// real number. This will be the standard error of the prior distribution of the standard deviation of the error terms  
    real<lower=0> scale_diag_D;// real number. This will be the standard error of the prior distribution of the standard deviation of the random terms 
    real<lower=0> lkj_shape;//real number. This will be the scale of the cholesky distribution
}
 
parameters {
    vector[ncx1] betas1;// coefficients of the variables in the fixed part 
    real<lower = 0> sigma1;// standard deviation of the error terms
    matrix[n, n_RE] b;// random terms b 
    vector<lower = 0>[n_RE] L_var_D;// vector with n_RE elements. It will contain the variances of the random effects b
    cholesky_factor_corr[n_RE] L_corr_D;// cholesky factor correlation matrix
}
 
transformed parameters {
    vector[N1] y1_hat = X1 * betas1; // fixed part 
    for (j1 in 1:N1) {
        for (q1 in 1:2) {
            y1_hat[j1] += Z_1[j1,q1] * b[id1[j1], RE_ind1[q1]];//estimated values 
        }
    }
}

  // lets define now the prior distributions 
  
model {
    matrix[n_RE, n_RE] L_D;
    vector[n_RE] mu;
    L_D = diag_pre_multiply(L_var_D, L_corr_D);
    L_var_D ~ cauchy(0, scale_diag_D);
    L_corr_D ~ lkj_corr_cholesky(lkj_shape);
    mu = rep_vector(0, n_RE);
    for (i in 1:n) {
        b[i, ] ~ multi_normal_cholesky(mu, L_D);
    }
    for (k1 in 1:ncx1) {
        betas1[k1] ~ normal(0.0, scale_betas1);
    }
    sigma1 ~ cauchy(0, scale_sigmas);
    y1 ~ normal(y1_hat, sigma1);
}
 
generated quantities {
    matrix[n_RE, n_RE] D;
    D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';// cholesky decomposition. We calculate the 
    //variance-covariance matrix of the random terms b 
}
