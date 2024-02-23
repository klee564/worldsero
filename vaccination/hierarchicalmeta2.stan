data {
  int<lower=0> J1;
  int<lower=0> J2;
  real y1[J1];  // num cases, treatment
  real<lower=0> v1[J1];  // num successes, treatment
  real y2[J2];  // num cases, treatment
  real<lower=0> v2[J2];  // num successes, treatment
  int<lower = 1> K;               // Number of groups.
  int<lower = 1, upper = K> g1[J1]; // Vector of group assignments.
  int<lower = 1, upper = K> g2[J2]; // Vector of group assignments.
}

parameters {
  //real theta[J];      // per-trial treatment effect
  real mu[K];            // mean treatment effect
  real<lower=0> dmu[K];            // mean treatment effect
  real<lower=0> tau;  // deviation of treatment effects
  real mu0;
  real<lower=0> dmu0;
  real<lower=0> kappa;  // deviation of treatment effects
  real<lower=0> pkappa;  // deviation of treatment effects
}

transformed parameters{

  real pmu[K];
  for(i in 1:K){
    pmu[i] = mu[i]+dmu[i];
  }

}
model {
  real gamma1[J1];
  real gamma2[J2];
  for(k in 1:K){
    mu[k] ~ normal(mu0,kappa);
    pmu[k] ~ normal(mu0+dmu0,pkappa);
  }
  for(j in 1:J1){
    y1[j] ~ normal(pmu[g1[j]],sqrt(tau^2+v1[j]));
    gamma1[j] = (tau/(v1[j]^2+tau^2))^2;
  }
  for(j in 1:J2){
    y2[j] ~ normal(mu[g2[j]],sqrt(tau^2+v2[j]));
    gamma2[j] = (tau/(v2[j]^2+tau^2))^2;
  }
  //mu0 ~ normal(0, 10);
  //dmu0 ~ normal(0, 10);
  //tau ~ cauchy(0, 5);
  target += log( sum(gamma1 ) + sum(gamma2))/2;
  //kappa ~ cauchy(0, 5);
  //pkappa ~ cauchy(0, 5);
  //y ~ normal(theta, sigma);
  //theta ~ normal(mu, tau);

}
