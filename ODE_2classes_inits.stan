// two age classes, with inital values that are drawn from exponential distribution, allowing for scaling of the begin size. [12.12]
functions {
    real[] SIR(real t,real[] y,real[] params,real[] x_r,int[] x_i) {
    
    real dydt[7]; // S, I and R for 3 groups, then C (cumulative incidence)
    real ntemp; // total population
    real nI; // total infectious
    real gamma;
    real beta[2];

    gamma = 1/(params[3]/7);
    beta[1] = params[1]*7;
    beta[2] = params[2]*7;

    ntemp = sum(y[1:6]);
    nI = y[2]+y[5];
    
    dydt[1] = -beta[1]  * y[1] * nI/ntemp; 
    dydt[2] = beta[1]  * y[1] * nI/ntemp - gamma * y[2];
    dydt[3] = gamma * y[2];

    dydt[4] = -beta[2] * y[4] * nI/ntemp; 
    dydt[5] = beta[2]   * y[4] * nI/ntemp - gamma * y[5];
    dydt[6] = gamma * y[5];

    dydt[7] = beta[1] * y[1] * nI/ntemp + beta[2]  * y[4] * nI/ntemp; 

    return dydt;
    }
  }
    
  data {
    int agegroups;
    int positives[agegroups];
    int N[agegroups];
    int inference;
    
    int<lower = 1> n_obs; // Number of weeks sampled
    int<lower = 1> n_params; // Number of model parameters
    int<lower = 1> n_difeq; // Number of differential equations in the system
    int y[n_obs]; // weekly observed incidence
    int tmax;
    real t0; // starting for y0[1], y0[2], y0[3] (S,I,R) point (zero)
    real ts[tmax];
    real initialsize[agegroups]; // initial population size at time zero
    //prior values
	real prior_beta1; // rate (days^-1)
    real prior_beta2; // rate (days^-1)
    real prior_gamma[2]; // rate (days^-1)
    real prior_rho[2]; // proportion
	real prior_init_inf; // exp
  }
    
  transformed data {
    real x_r[0];
    int x_i[0];
  }
    
  parameters {
    real<lower = 0> params[n_params]; // Model parameters
    real<lower = 0,upper=1> rho;
    real<lower = 0> sigma;
	real<lower = 0>init_inf;
  }
  transformed parameters{
    real y_hat[tmax, n_difeq]; // Output from the ODE solver
    real C[n_obs];
    real prop_end[agegroups]; //compute end proportions
    real y0[7];

    y0[1] = initialsize[1]-init_inf/2;
    y0[2] = init_inf/2;
    y0[3] = 0;

    y0[4] = initialsize[2]-init_inf/2;
    y0[5] = init_inf/2;
    y0[6] = 0;
    y0[7] = 0;
	
    y_hat = integrate_ode_rk45(SIR, y0, t0, ts, params, x_r, x_i, 1.0E-6, 1.0E-6, 1.0E3); //set max depth 
    // JULIEN: changed settings according to https://discourse.mc-stan.org/t/speeding-up-a-stochastic-epidemic-model-with-odes/941/8
    prop_end[1] = y_hat[tmax,3]/(y_hat[tmax,1]+y_hat[tmax,2]+y_hat[tmax,3]);
    prop_end[2] = y_hat[tmax,6]/(y_hat[tmax,4]+y_hat[tmax,5]+y_hat[tmax,6]);

    for (n in 1:n_obs)
      C[n]=rho*(y_hat[n,7]- (n==1 ? 0 : (y_hat[n,7]>y_hat[n-1,7] ? y_hat[n-1,7] : 0) )); //if n=1, substract 0, else 
      // JULIEN: found the problem with negative values: if the epidemic dies out (if a low beta is sampled for instance)
      // then C is stable, and with some numerical errors the difference between successive C can become negative
  }
    
  model {
    params[1]~exponential(prior_beta1); 
    params[2]~exponential(prior_beta2);
    params[3]~gamma(prior_gamma[1],prior_gamma[2]); 
	
    rho~beta(prior_rho[1],prior_rho[2]);
    sigma~cauchy(0,3);
	init_inf~exponential(prior_init_inf); // initial amount of infected people (divided over I_1, I_2)    
    
    if (inference==1){
      //y~poisson(C);
      sqrt(y)~normal(sqrt(C), sigma);
      positives[1]~binomial(N[1], prop_end[1]);
      positives[2]~binomial(N[2], prop_end[2]);

    }
  }
  generated quantities{
    real pred_Y[n_obs];
    real pred_p[agegroups];
	vector[n_obs+agegroups] log_lik;
 
    for (n in 1:n_obs)
      pred_Y[n]=(normal_rng(sqrt(C[n]), sigma))^2;
    for (n in 1:agegroups)
      pred_p[n] = binomial_rng(N[n], prop_end[n]);
	  
	 for (n in 1:n_obs) {
		log_lik[n] = normal_lpdf(sqrt(y[n]) | sqrt(C[n]), sigma);
		}
	 for (n in 1:agegroups)
      log_lik[n_obs+n] = binomial_lpmf(positives[n] | N[n], prop_end[n]);
  }
