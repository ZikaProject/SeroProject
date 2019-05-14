/* empty stan model that generates quantities.
 first attempt at ABM. */
functions {
  real RateToProb01(real rate){ //returns probability per week based on daily rate.
	return(1-exp(-(rate*0.1)));
  }
 // real RateToProbwk(real rate){ //returns probability per week based on daily rate.
//	return(1-exp(-(rate*7)));
 // }
  int test_func();
  // death function: returns probability of dying that week for specific age.
  real death(real age, vector mort_rates){ // return probability of dying for an individual per age
    {
      real m_rate;
      if(age<1){ // probability of death <1 mort_rate[1].
        m_rate = mort_rates[1]/365.25; // yearly rate to daily
      }else{ // use polynomial function that was fit to data.
        m_rate = exp(mort_rates[2]+(mort_rates[3]*age)+(mort_rates[4]*age^2)+(mort_rates[5]*age^3))/365.25;
      }
      return(RateToProb01(m_rate)); 
    }
  }
  
  // birth function: returns probability of birth that week for specific age.
  real birth(real age, real b_rate){
    if(age>=15 && age<50){
      return(RateToProb01(b_rate/(70*365.25))); //birth rate per lifetime (reproductive timespan assume 35y)
    }else{
      return(0);
    }
  }
   
  //return probability of infection for that week for specific age.
  real infect(real age, real beta1, real beta2, real num_inf, real num_tot){
    if(age<15){
      return(RateToProb01(beta1*(num_inf/num_tot))); //beta: rate per day
    }else{
      return(RateToProb01(beta2*(num_inf/num_tot)));
    }
  }
 
 //return probability of infection for that week for specific age (now still age-indep)
   real recover(real age, real gamma){
    return(RateToProb01(gamma)); // gamma is in days
  }
  
  //return probability of loss of immunity:
  // todo
  
  // classify age:
  vector ageclass(real age){
    vector[4] age_return;
    if(age<15){
      age_return=to_vector({1,0,0,0}); //children
    }else if(age>=15 && age <30){ // young, childbearing
      age_return=to_vector({0,1,0,0}); 
    }else if(age>=30 && age <50){ // old, childbearing
      age_return=to_vector({0,0,1,0});
    }else{
      age_return=to_vector({0,0,0,1}); // 50 and up
    }
    return(age_return);
  }
  
  //probability of S, I, R for initial population
  real initialize_SIR(real age, real prop1, real prop2){ //S status inital population
	if(age<15){
		return(prop1); //observed proportion that were infected (R) from data.
	}else{
		return(prop2);
	}
  } 

}
data {
  int N0; // initial number of individuals in 
  int nend; // maximum number of people allowed in 
  int units; // number of simulated time units
  int runN; // nth simulation.
  int introduction_wk; //wk of introduction of infection
  int introduction_n; //n introduced
  real age_init[nend]; // initial ages of the N0 individuals
  real birthrate; // the birth rate (number of neonates per women). 
  vector[5] mort_rates; // [1] mort_rate <1, [2:5] coefficients of polynomial describing >
  real beta1[4000]; //https://discourse.mc-stan.org/t/is-it-possible-to-access-the-iteration-step-number-inside-a-stan-program/1871/14
  real beta2[4000];
  real gamma[4000];
  real prop1[4000];
  real prop2[4000];
}
  parameters {
   real<lower = 0> params[3]; 
    
}
transformed parameters {
}


model {

}

generated quantities{

vector[nend] ind_alive;
vector[nend] ind_age;
vector[nend] ind_S;
vector[nend] ind_I;
vector[nend] ind_R;
vector[nend] ind_C;
//vector[nend] agegrp[4];
vector[nend] agegrp1;
vector[nend] agegrp2;
vector[nend] agegrp3;
vector[nend] agegrp4;

//vector[4] ageret;

real alive[4174]; // calculate number of people alive for each week.
real infected[4174];
real population; // if population=1: ageing, death, birth, otherwise switched off.

//real n_alive[units+1];
//real n_infected[units+1];
real recovered[4174];
real susceptible[4174];
real cumind[4174];
real recovered_age1[4175];
real recovered_age2[4175];

real c_inf; //cumulative number of new infected people.
real c_inf1; // per age group
real c_inf2;
real c_inf3;
real c_inf4;

//keep track of population at start of the outbreak (at time introduction_wk)
real alive_start;

real S_start;
real S_start1;
real S_start2;
real S_start3;
real S_start4;

real R_start;
real R_start1;
real R_start2;
real R_start3;
real R_start4;
real alive_n;
real infected_n;

real prop_c_inf;
int total_N;
int infectind;
int b;
int prev;
int iteration;
int iteration_n;
int cnt_wk;

population = 1;
cnt_wk = 1;
total_N = N0; // initial population size from data.
iteration = test_func();
iteration_n = iteration + (runN*10); // runN is nth run of 10 iterations.
// initialise age groups, SIR and alive status:
//print("initialize:");
for (i in 1:nend){
	ind_age[i] = age_init[i];
	ind_C[i] = 0;
	alive_n = N0;
	infected_n = 0;
   if(ind_age[i]<15){
      	agegrp1[i] = 1;
		agegrp2[i] = 0;
		agegrp3[i] = 0;
		agegrp4[i] = 0; 
    }else if(ind_age[i]>=15 && ind_age[i] <30){ // young, childbearing
     	agegrp1[i] = 0;
		agegrp2[i] = 1;
		agegrp3[i] = 0;
		agegrp4[i] = 0; 
    }else if(ind_age[i]>=30 && ind_age[i] <50){ // old, childbearing
     	agegrp1[i] = 0;
		agegrp2[i] = 0;
		agegrp3[i] = 1;
		agegrp4[i] = 0; 
    }else{
     	agegrp1[i] = 0;
		agegrp2[i] = 0;
		agegrp3[i] = 0;
		agegrp4[i] = 1; 
    }
	// SR initialisation based on N and positive (per draw from probability)?
	if(i<=N0){ //active people
		ind_R[i] = bernoulli_rng(initialize_SIR(ind_age[i], prop1[iteration_n], prop2[iteration_n]));
		ind_S[i] = 1-ind_R[i];
		ind_I[i] = 0;
		ind_alive[i] = 1;
	}else{ //spare people
		ind_S[i] = 1;
		ind_R[i] = 0;
		ind_I[i] = 0;
		ind_alive[i] = 0;
	}		
}

//calculate initial population summaries
alive[1] = sum(ind_alive);
infected[1] = sum(ind_alive .* ind_I);
cumind[1]=0;
recovered[1] = sum(ind_alive .* ind_R);
//prop_recovered[1]= sum(ind_alive .* ind_R)/sum(ind_alive);
recovered_age1[1] = sum(ind_alive .* ind_R .*agegrp1);
recovered_age2[1] = sum(ind_alive .* ind_R .*agegrp2);


//c_inf = 0;
prop_c_inf = 0;

//prop_recovered_age1[1] = 0;
//mean_age[1] = sum(ind_alive .* ind_age)/sum(ind_alive);
// simulate for units of time:
   //print("THIS IS ITERATION:", iteration, " beta1:", beta1[iteration], " prop1:", prop1[iteration]);
   for (n in 1:units){    
	for (i in 1:total_N){ // per week go through all individuals		
		// check if this person dies this week
		if(ind_alive[i]==1){ // if alive:			  
			  
			
				// birth:
				if(bernoulli_rng(birth(ind_age[i],birthrate))==1){
					total_N=total_N+1;
					// 'activate' newborn:
					ind_alive[total_N]=1; // make individual total_N alive.
				}
				// aging:
				ind_age[i] = ind_age[i]+(1/(365.25*10)); 
				//calculate age-groups, can probably be coded more efficiently.
			   if(ind_age[i]<15){
					agegrp1[i] = 1;
					agegrp2[i] = 0;
					agegrp3[i] = 0;
					agegrp4[i] = 0; 
				}else if(ind_age[i]>=15 && ind_age[i] <30){ // young, childbearing
					agegrp1[i] = 0;
					agegrp2[i] = 1;
					agegrp3[i] = 0;
					agegrp4[i] = 0; 
				}else if(ind_age[i]>=30 && ind_age[i] <50){ // old, childbearing
					agegrp1[i] = 0;
					agegrp2[i] = 0;
					agegrp3[i] = 1;
					agegrp4[i] = 0; 
				}else{
					agegrp1[i] = 0;
					agegrp2[i] = 0;
					agegrp3[i] = 0;
					agegrp4[i] = 1; 
				}
				// death:
				ind_alive[i] = 1-bernoulli_rng(death(ind_age[i], mort_rates));
			
			  // recovery:
			//infectious actions
			  if(ind_I[i]==1 && bernoulli_rng(recover(ind_age[i], (1/gamma[iteration_n])))==1){
				  ind_I[i]=0;
				  ind_R[i]=1;	
				  //print("individual ", i , " recovered in wk ",n,", alive: ", alive[n]," & infected:", infected[n])
			  }			  
			  // infection, probability of infection depends on age, beta1,2 number of infected and total alive:
			  if(ind_S[i]==1 && bernoulli_rng(infect(ind_age[i], beta1[iteration_n], beta2[iteration_n], infected_n, alive_n))==1){
			    // add if inf>0 not
				  ind_S[i]=0;
				  ind_I[i]=1;
				  //c_inf = c_inf+1;
				ind_C[i] = 1;
				  //print("individual ", i , " got infected in wk ",n,", alive: ", alive[n]," & infected:", infected[n])
			  } // still have births in 
			  
			


			}
			//ind[i,2] = death(ind[i,2]);			
	}	
	// code below introduces [introduction_n] infections at [introduction_wk], coule be wrapped in a function.	
	
	 if(n>1){
		prev=n-1;
		// whole year:
		if(floor(prev/3652)!=floor(n/3652)){
		  b=0;
		  print("next year: ", n/3652);
		  while (b < introduction_n) { // introduce [introduction_n] new infections among alives
			  infectind = categorical_rng(rep_vector(inv(total_N), total_N)); //draw random person from total_N
			  if(ind_alive[infectind]==1 && ind_S[infectind]==1){ // infect when alive and susceptible, and move to the next
				ind_S[infectind]=0;
				ind_I[infectind]=1;
				ind_C[infectind]=1; // introductions are also counted as incidenct cases.
				b=b+1;
				
			  }
				  
			}
			print("INTRODUCTION", sum(ind_alive), "=alive ", sum(ind_C), "=C", iteration_n, "=iteration");
		}
		// for 4172 weeks: 
		//4175
		if(floor(prev/70)!=floor(n/70)){
		  //print("next wk: ", n/70, " equals: ", cnt_wk);
		  
		  
		  //properties to keep per week:
		  infected[cnt_wk]=sum(ind_C);//sum(ind_alive .* ind_C);
		  if(cnt_wk>1){
			cumind[cnt_wk]=infected[cnt_wk]-infected[cnt_wk-1];
			}
		  susceptible[cnt_wk]=sum(ind_alive .* ind_S);
		  recovered[cnt_wk]=sum(ind_alive .* ind_R);
		  cnt_wk=cnt_wk+1;
		  
		}
	}
 // at time [introduction_wk] wks

	   /*alive_start = sum(ind_alive);
		S_start = sum(ind_alive .* ind_S);
		S_start1 = sum(ind_alive .* ind_S .*agegrp1);
		S_start2 = sum(ind_alive .* ind_S .*agegrp2);
		S_start3 = sum(ind_alive .* ind_S .*agegrp3);
		S_start4 = sum(ind_alive .* ind_S .*agegrp4);

		R_start = sum(ind_alive .* ind_R);
		R_start1 = sum(ind_alive .* ind_R .*agegrp1);
		R_start2 = sum(ind_alive .* ind_R .*agegrp2);
		R_start3 = sum(ind_alive .* ind_R .*agegrp3);
		R_start4 = sum(ind_alive .* ind_R .*agegrp4);
*/
	  // calculate S, R at start of epidemic.

    //print("alive: ", alive[n]);
	
	
	/*
	for (i in 1:nend){
    	ind_alive[i] = bernoulli_rng(probabs[i]) .* ind_alive[i];
		  // keep vector of probabilities in local environment in function?
		  // have function return vector of length nend with 0's or 1's.
		  // move bernoulli_rng to function
		}*/	
	// calculate totals
	
	alive_n = sum(ind_alive);
	infected_n = sum(ind_alive .* ind_I);
	//c_tot = c_tot + (infected[n+1] - infected[n])
	
	c_inf = sum(ind_C);
	c_inf1 = sum(ind_C .* agegrp1);
	c_inf2 = sum(ind_C .* agegrp2);
	c_inf3 = sum(ind_C .* agegrp3);
	c_inf4 = sum(ind_C .* agegrp4);
	//prop_c_inf = c_inf/alive[n+1];
	// c_inf per age group will be informative to see where the outbreak proportionally hits hardest (likely young child bearing age, group2)
	
	//recovered[n+1] = sum(ind_alive .* ind_R);
	//recovered_age1[n+1] = sum(ind_alive .* ind_R .*agegrp1);
	//recovered_age2[n+1] = sum(ind_alive .* ind_R .*agegrp2);
	//prop_recovered[n+1]= sum(ind_alive .* ind_R)/sum(ind_alive);
	// per age group:
	//recovered_age1[n+1] = sum(ind_alive .* ind_R .* agegrp[1,]);
	//prop_recovered_age1[n] = sum(ind_alive .* ind_R .* agegrp[1,])/sum(ind_alive .* agegrp[1,]);
	//prop_recovered_age2[n] = sum(ind_alive .* ind_R .* agegrp[2,])/sum(ind_alive .* agegrp[2,]);
	//relprop_recovered_age2[n] = sum(ind_alive .* ind_R .* agegrp[2,])/(sum(ind_alive .* ind_R .* agegrp[1,])+sum(ind_alive .* ind_R .* agegrp[2,])+sum(ind_alive .* ind_R .* agegrp[3,])+sum(ind_alive .* ind_R .* agegrp[4,]));
	/*recovered_age2[n+1] = 
	recovered_age3[n+1] = 
	recovered_age4[n+1] = 
	*/
	
	//mean_age[n+1] = sum(ind_alive .* ind_age)/sum(ind_alive);
   }
   
}

 
