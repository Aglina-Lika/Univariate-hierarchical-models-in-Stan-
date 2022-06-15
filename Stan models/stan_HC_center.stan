data {
    int n;//total number of subjects 
    int n_RE; // total numer of random effects
    int N1;  // total number of observations 
    int ncx1; // total number of variables in the fixed part
    int id1[N1]; // vector of repeated ids
    int RE_ind1[2]; // vector with 2 elements
    vector[N1] y1; // the dependent variable 
    matrix[N1, 2] Z_1; // the design matrix of the random effects when we have raw variables 
    matrix[n, ncx1] Xhc1; // the design matrix of the fixed part when we have raw variables 
    int ZrowsStart1[n];  // vector of n elements. This vector will have the position where is the first observation of each patient 
    int ZrowsEnd1[n]; // vector of n elements. This vector will have the position where is the last observation of each patient 
    matrix[2, N1] Zv1; // transpose design matrix of the random effects when we have raw variables 
    real<lower=0> scale_betas1; // real number. This will be the standard error of the prior distribution of the coefficients 
    real<lower=0> scale_sigmas; // real number. This will be the standard error of the prior distribution of the standard deviation of the error terms  
    real<lower=0> scale_diag_D; // real number. This will be the standard error of the prior distribution of the standard deviation of the random terms 
    real<lower=0> lkj_shape; //real number. This will be the scale of the cholesky distribution
    matrix[n, 4] XhcC1; // the design matrix of the fixed part when we have centred variables 
    vector[2] means_Xhc1; // means of the raw variables in the fixed part 
    matrix[N1,2] Zc1; // the design matrix of the random part when we have centred variables
    real means_Z1; // means of the raw variables in the random part 
    matrix[N1,2] Ztinv1; // the transposed inversed matrix of the random part when we have raw variables 
    matrix[2,N1] Zinv1;// the inversed matrix of the random part when we have raw variables 
}
 
parameters {
    vector[4] temp_betas1; // coefficients of the variables in the fixed part 
    real<lower = 0> sigma1; // standard deviation of the error terms 
    matrix[n, n_RE] u; // random terms u 
    vector<lower = 0>[n_RE] L_var_D; // vector with n_RE elements. It will contain the variances of the random effects u
    cholesky_factor_corr[n_RE] L_corr_D; // cholesky factor correlation matrix
}
 
transformed parameters {
    vector[N1] eta1; // it will include the expected values of the model 
    matrix[n, n_RE] mu_u;
    for (i in 1:n) {
        mu_u[i, 1] = XhcC1[i, 1] * temp_betas1[1] + XhcC1[i, 3] * temp_betas1[3] + XhcC1[i, 4] * temp_betas1[4]; // it will give the mean
        // of the multivariate normal distribution of the random intercept of the u 
        mu_u[i, 2] = XhcC1[i, 2] * temp_betas1[2]; // it will give the mean
        // of the multivariate normal distribution of the random slope of the u
    }
    for (j1 in 1:N1) {
        eta1[j1] = Zc1[j1, 1] * u[id1[j1], 1] + Zc1[j1, 2] * u[id1[j1], 2]; // estimated values 
    }
}

 // lets define now the prior distributions 
 
model {
    matrix[n_RE, n_RE] L_D;
    L_D = diag_pre_multiply(L_var_D, L_corr_D);
    L_var_D ~ cauchy(0, scale_diag_D); // we define a cauchy distribution 
    L_corr_D ~ lkj_corr_cholesky(lkj_shape);
    for (i in 1:n) {
        u[i, ] ~ multi_normal_cholesky(mu_u[i, ], L_D);
    }
    for (k1 in 2:4) {
        temp_betas1[k1] ~ normal(0.0, scale_betas1);
    }
   sigma1 ~ cauchy(0, scale_sigmas);
    y1 ~ normal(eta1, sigma1);
}
 
generated quantities {
    vector[ncx1] betas1; // these are the transformed coefficients 
    real Intercept1 = temp_betas1[1] - means_Xhc1[1] * temp_betas1[3] - means_Xhc1[2] * temp_betas1[4] - means_Z1 * temp_betas1[2];// we tranform back the intercept 
    matrix[n_RE, n_RE] D; // variance covariance matrix of the random effects b
    matrix[n_RE, n_RE] tmp_D; // variance covariance matrix of the random effects u
    matrix[n, n_RE] b; // these will be the random effects b 
    matrix[n, n_RE] tmp_b; // we define a matrix of dimension nxn_RE. This will have in the columns the random intercept and random slope 
    // of the random terms b 
    betas1= append_row(Intercept1, temp_betas1[2: 4]); // we define the vector of coefficient that the STAN code will return
    tmp_D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)'; // cholesky decomposition. We calculate the 
    //variance-covariance matrix of the random terms u 
    tmp_b = u - mu_u; // we calcuate the random terms b
    D =  Zinv1 * Zc1 * tmp_D * (Zc1)' * Ztinv1; # we transform back the variance covariance matrix of the random terms 

        b[, 1] = tmp_b[, 1]-tmp_b[, 2]*means_Z1; # we trasform back the random intercept 

        b[,2]=tmp_b[,2]; # no need to transform back the random slope 
}
