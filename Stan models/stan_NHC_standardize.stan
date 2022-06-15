data {
    int n;//total number of subjects 
    int n_RE;// total numer of random effects
    int N1;// total number of observations
    int ncx1;// total number of variables in the fixed part
    int id1[N1];// vector of repeated ids
    int RE_ind1[2];// vector with 2 elements
    vector[N1] y1;// the dependent variable 
    matrix[N1, 2] Zs1;// the design matrix of the random part when we have standardized variables 
    real<lower=0> scale_betas1;// real number. This will be the standard error of the prior distribution of the coefficients 
    real<lower=0> scale_sigmas;// real number. This will be the standard error of the prior distribution of the standard deviation of the error terms  
    real<lower=0> scale_diag_D;// real number. This will be the standard error of the prior distribution of the standard deviation of the random terms 
    real<lower=0> lkj_shape;//real number. This will be the scale of the cholesky distribution
    int ZrowsStart1[n];// vector of n elements. This vector will have the position where is the first observation of each patient 
    int ZrowsEnd1[n];// vector of n elements. This vector will have the position where is the last observation of each patient 
    matrix[2, N1] Zv1;// the trasposed design matrix of the random part when we have raw variables 
    real SDs_Z1;// the vector contains the standard deviations of the variables in the random part 
    real mean_sd_Z1; //the vector contains the means/sds of the variables in the random part 
    matrix[N1, 4] Xs1; // the design matrix of the fixed part when we have standardized variables 
    vector[3] SDs_X1;// the vector contains the standard deviations of the variables in the fixed part 
 
    vector[3] mean_sd_X1;//the vector contains the means/sds of the variables in the fixed part 
    matrix[N1,2] Ztinv1;//// the trasposed inverse design matrix of the random part when we have raw variables 
    matrix[2,N1] Zinv1;// the inverse design matrix of the random part when we have raw variables 
}
 
parameters {
    vector[4] temp_betas1;// coefficients of the variables in the fixed part 
    real<lower = 0> sigma1;// standard deviation of the error terms
    matrix[n, n_RE] v;// random terms v 
    vector<lower = 0>[n_RE] L_var_D;// vector with n_RE elements. It will contain the variances of the random effects v
    cholesky_factor_corr[n_RE] L_corr_D;// cholesky factor correlation matrix
}
 
transformed parameters {
    vector[N1] y1_hat = Xs1 * temp_betas1;// fixed part
    for (j1 in 1:N1) {
        for (q1 in 1:2) {
            y1_hat[j1] += Zs1[j1,q1] * v[id1[j1], RE_ind1[q1]];//estimated value
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
        v[i, ] ~ multi_normal_cholesky(mu, L_D);
    }
    for (k1 in 2:4) {
        temp_betas1[k1] ~ normal(0.0, scale_betas1);
    }
    sigma1 ~ cauchy(0, scale_sigmas);
    y1 ~ normal(y1_hat, sigma1);
}
 
generated quantities {
    vector[ncx1] betas1;// it will contain the transformed intercept and the 
      // other coefficients that do not need to be trasformed back 
    real Intercept1 = temp_betas1[1] - dot_product(mean_sd_X1, temp_betas1[2:4]);// we transform 
  // back the intercept 
    vector[3] betas1_part;// it will contain the transformed slopes (coefficients)
    matrix[n_RE, n_RE] D;// the transformed variance covariance matrix  
    matrix[n, n_RE] b; // random effects b 
        matrix[n_RE, n_RE] tmp_D;//variance covariance matrix of the random effects v
//with the for loop we transform back also, the other coefficients of the fixed part 
    for (i in 1:3) {
        betas1_part[i] = temp_betas1[i+1] / SDs_X1[i];
    }
    betas1= append_row(Intercept1, betas1_part);


    tmp_D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';//
    D =  Zinv1 * Zs1 * tmp_D * (Zs1)' * Ztinv1;//

 
        b[, 1] = v[, 1]-v[, 2]*mean_sd_Z1; // we transform back the random intercept
  
        b[,2]=v[,2]/SDs_Z1; // we transform back the random slope 
}
