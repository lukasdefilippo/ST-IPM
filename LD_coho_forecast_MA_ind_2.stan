data{
	//Number of years (total, includes several missing years for some stocks)
	int <lower=0> n_year;

	//Number of total populations
	int <lower=0> n_pop;

	//Number of populations with return data
	int <lower=0> n_pop_tot;

	//Which populations possess return data
	int <lower=0> pop_tot[n_pop_tot];

	//Length of the return data vectors
	int <lower=0> n_tot;

	//Vectors of all return data across all populations
	real <lower=0> tot_dat [n_tot];

	//Vectors of the indices identifying which years are those with non-NA data for the return data
	int <lower=0> tot_true [n_tot];

	//Paired vectors of slice points indicating the beginning, and end of the data for a particular population
	int <lower=0> slice_tot_start [n_pop_tot];
	int <lower=0> slice_tot_end [n_pop_tot];

	
}parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  	//Initial adult abundance
  	vector <lower=0> [n_pop] sigma_obs;

	// Population-specific process error variance
	vector <lower=0> [n_pop] sigma_proc;

	vector <lower=-0.99999, upper = 0.99999> [n_pop] phi;

	matrix [n_year,n_pop] proc_dev;

	// log of initial adult abundance
  	vector [n_pop] init;

  	vector [n_pop] mu;

}transformed parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// State estimates of log adult returns

	// Populaiton-specific process errors
	matrix [n_year,n_pop] log_adult_est;

    // Populaiton-specific process errors
	matrix [n_year, n_pop] epsilon;

	matrix [n_year, n_pop] pred;

  	// State estimates of adult returns
  	matrix <lower=0> [n_year,n_pop] adult_est;

	//Population-specific random walk of adult returns
	for(i in 1:n_pop){
		for(y in 1:n_year){
			if(y==1){
				pred[y,i] = init[i];

				log_adult_est[y,i] = pred[y,i] + proc_dev[y,i];

				epsilon[y,i] = log_adult_est[y,i] - pred[y,i];

			}else{
				
				pred[y,i] = mu[i] + phi[i] * epsilon[y - 1, i];

				log_adult_est[y,i] = pred[y,i] + proc_dev[y, i];

				epsilon[y,i] = (log_adult_est[y,i] - mu[i] - phi[i] * epsilon[y - 1,i]);

			}
		}
	}

	//exponentiate adult returns from estimation in log space
	adult_est = exp(log_adult_est);

}model {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	
 	//Process deviations drawn from population-specific error distributions
 	for(i in 1:n_pop){
 		proc_dev[,i] ~ normal(0, sigma_proc[i]);
 	}

 	//Vague prior for initial log abundance
 	mu ~ normal(0,20);

 	sigma_obs ~ normal(0,5);

 	init ~ normal(0,20);

 	// Prior for 1st year adult return variation
 	sigma_proc~ normal(0,5);

	//Likelihood
	for(i in 1:n_pop_tot){
		tot_dat[slice_tot_start[i]:slice_tot_end[i]] ~ lognormal(log_adult_est[tot_true[slice_tot_start[i]:slice_tot_end[i]], pop_tot[i]], sigma_obs[i]);
	}

} generated quantities {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Vector of adult return forecasts
    vector[n_pop] adult_pred;

    for(i in 1:n_pop){
		adult_pred[i] = exp(log_adult_est[n_year, i]);
	}

}
