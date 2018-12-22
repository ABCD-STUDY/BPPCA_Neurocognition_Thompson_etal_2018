data{  int<lower=1> N; 	// number of subjects in dataset  int<lower=1> P; 	// dimension of data
  int<lower=1> D; 	// dimension of latent space  int<lower=1> ns;	// number of sites
  int<lower=1> nf;	// number of families with more than one member
  int<lower=1> Site[N]; // site ids
  int<lower=0> Fam[N];	// family ids
  vector[P] Y[N]; 	// responses
}
transformed data{
  vector[P] eta;
  eta = rep_vector(0,P);
}
parameters{  vector[D] Theta[N];  matrix[P,D] Lambda;
  matrix[D,ns] A_site;
  matrix[D,nf] B_fam;  matrix[P,ns] C_site;
  matrix[P,nf] D_fam;  real<lower=0> sigma2_a;
  real<lower=0> sigma2_b;
  real<lower=0> sigma_c;
  real<lower=0> sigma_d;
  real<lower=0> sigma_eps;}
transformed parameters{
  vector[D] mu[N];
  vector[P] nu[N];
  for(n in 1:N){
    if(Fam[n] == 0) mu[n] = col(A_site,Site[n]);    if(Fam[n]>0) mu[n] = col(A_site,Site[n]) + col(B_fam,Fam[n]);
  }
  for(n in 1:N){
    if(Fam[n] == 0) nu[n] = col(C_site,Site[n]);    if(Fam[n]>0) nu[n] = col(C_site,Site[n]) + col(D_fam,Fam[n]);
  }
}model{  sigma2_a ~ uniform(0,1);
  sigma2_b ~ uniform(0,1-sigma2_a);
  sigma_c ~ cauchy(0,1);
  sigma_d ~ cauchy(0,1);
  sigma_eps ~ cauchy(sigma_c^2+sigma_d^2,1);
  to_vector(A_site) ~ normal(0,sqrt(sigma2_a));
  to_vector(B_fam) ~ normal(0,sqrt(sigma2_b));
  to_vector(C_site) ~ normal(0,sigma_c);
  to_vector(D_fam) ~ normal(0,sigma_d);
  for(d in 1:D)
    Lambda[,d] ~ normal(0, 1);  for(n in 1:N){    if(Fam[n] == 0){
	Theta[n] ~ multi_normal(mu[n], diag_matrix(rep_vector(1-sigma2_a,D)));
    	Y[n] ~ multi_normal(Lambda * Theta[n] + nu[n], diag_matrix(rep_vector(sigma_eps^2-sigma_c^2,P)));
    }
    if(Fam[n] > 0){
	Theta[n] ~ multi_normal(mu[n], diag_matrix(rep_vector(1-sigma2_a-sigma2_b,D)));    	Y[n] ~ multi_normal(Lambda * Theta[n] + nu[n], diag_matrix(rep_vector(sigma_eps^2-sigma_c^2-sigma_d^2,P)));
    }  }}
generated quantities{
  //matrix[P,P] Eig_vec_lambda;
  //vector[P] Eig_val_lambda;
  //matrix[P,P] W;
  cov_matrix[P] Q;
  vector[N] log_lik_marg;
  vector[N] log_lik;
  //Eig_vec_lambda = eigenvectors_sym(Lambda*Lambda');
  //Eig_val_lambda = eigenvalues_sym(Lambda*Lambda');
  Q = Lambda * Lambda' + diag_matrix(rep_vector(sigma_eps^2,P));
  //W = Lambda * Lambda';
  for (n in 1:N){
    log_lik_marg[n] = multi_normal_lpdf(Y[n] | eta, Q);
    if(Fam[n] == 0){
     	log_lik[n] = multi_normal_lpdf(Y[n] | Lambda * Theta[n] + nu[n], diag_matrix(rep_vector(sigma_eps^2-sigma_c^2,P)));
    }
    if(Fam[n] > 0){
     	log_lik[n] = multi_normal_lpdf(Y[n] | Lambda * Theta[n] + nu[n], diag_matrix(rep_vector(sigma_eps^2-sigma_c^2-sigma_d^2,P)));
    }  }
}
