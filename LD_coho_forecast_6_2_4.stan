data{
	//Number of years (total, includes several missing years for some stocks)
	int <lower=0> n_year;

	//Number of total populations
	int <lower=0> n_pop;

	//Number of populations of each data type
	int <lower=0> n_pop_smolt;
	int <lower=0> n_pop_esc;
	int <lower=0> n_pop_catch;
	int <lower=0> n_pop_MS;

	//Which populations possess have each data type
	int <lower=0> pop_smolt [n_pop_smolt];
	int <lower=0> pop_esc [n_pop_esc];
	int <lower=0> pop_catch [n_pop_catch];
	int <lower=0> pop_MS [n_pop_MS];

	//Number of populations whose marine survival data is from a nearby hatchery
	int <lower=0> n_hatchery;

	//Vector identifying which populations have marine survival estimates derived from a hatchery
	int <lower=0> hatchery[n_hatchery];
	int <lower=0> wild[n_pop_MS-n_hatchery];

	//Length of the smolt, harvest and escapement data vectors (the smolt, harvest, and escapement data for sll stocks are each read in as a single vector and sliced up among individual populations using a 
	//vector of cut points for each data type (see below) because stan does not accept ragged data structures)
	int <lower=0> n_smolt;
	int <lower=0> n_esc;
	int <lower=0> n_harvest;
	int <lower=0> n_MS;

	//Vectors of all smolt, harvest, and escapement data across all populations
	real <lower=0> smolt_dat [n_smolt];
	real <lower=0> esc_dat [n_esc];
	real <lower=0> harvest_dat [n_harvest];

	//Vectors of the number or recoveries as releases for the CWT marine survival data
	int MS_dat_x [n_MS];
	int MS_dat_N [n_MS];

	//Vectors of the indices identifying which years are those with non-NA data for the smolt, escapement, and harvest data
	int <lower=0> smolt_true [n_smolt];
	int <lower=0> esc_true [n_esc];
	int <lower=0> harvest_true [n_harvest];
	int <lower=0> MS_true [n_MS];

	//Paired vectors of slice points indicating the beginning, and end of the data for a particular population within the larger vector of that data type
	//The 'start' vector indicates the index when a particular population's data begins, and the end vector indicates when that population's data terminates
	// There is a start and end point for each population, so each vector is of length n_pop
	//Smolt data
	int <lower=0> slice_smolt_start [n_pop_smolt];
	int <lower=0> slice_smolt_end [n_pop_smolt];
	//Escapement data
	int <lower=0> slice_esc_start [n_pop_esc];
	int <lower=0> slice_esc_end [n_pop_esc];
	//Harvest data
	int <lower=0> slice_harvest_start [n_pop_catch];
	int <lower=0> slice_harvest_end [n_pop_catch];
	//Coded wire tag data
	int <lower=0> slice_MS_start [n_pop_MS];
	int <lower=0> slice_MS_end [n_pop_MS];

	//Observation error terms for the smolt, escapement, and harvest data (If fixed, supplied as data, otherwise estimated)
	//real <lower=0> sigma_smolt;
	real <lower=0> sigma_esc;
	//real <lower=0> sigma_catch;

	//An identity matrix used to multiply by the gamma matrix containing the means of the hyperdistributions for the S-R parameters to produce compatible dimsensions so that it can 
	//be added to the covariance matrix to implement the multivariate non-centered parameterization
	matrix[1, n_pop] u;

	//The eastings and northings for each population's marine entry point, expressed as an array of vectors with a dimension of 2 (easting + northings) and a length equal to the number of populations (n_pop)
  	matrix[n_pop,n_pop] dist;

  	//Stream distances for each habitat
  	vector <lower=0> [n_pop] stream_dist;

}transformed data {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Mean of the marine survival distribution (vector of 0s)
	vector[n_pop] mu_MS = rep_vector(0, n_pop);

}parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Uniform variable used to specify a re-parameterized cauchy prior on the concentration term (K) of the beta binomial likelihood for the CWT marine survival data
	real<lower=2.0, upper=500> k;

	//Mean and standard deviation for the hyperdistribution from which each individual population's Beverton-Holt error variance is drawn
	real mu_sigma_R;
	real <lower=0> sigma_sigma_R;

	//Scaling term for the non-centered parameterization of the hierarchically distributed recruitment error varainces
	vector[n_pop] z_sigma_R;

	//Matrix of annual recruitment deviations by population
    matrix [n_year,n_pop] r_dev;

    //Matrix of intitial spawning abundance for the first two recruitment events arising from unobserved spawning events, modelled as a latent parameter as in Fleishman et al. (2013)
    matrix [2, n_pop] log_adult_init;

	//Matrix of means for the shared distribution of the S-R parameters where row 1 is the mean of the productivity term (alpha) and row 2 is the mean of the asymptotic recruitment term (R_max)
	matrix[2, 1] gamma;

	//Cholesky factor of the correlation matrix for the S-R parameters
	cholesky_factor_corr[2] L_Omega;

	//Matrix of scaling factors for the multivariate implementation of the non-centered parameterization for the S-R parameters
	matrix[2, n_pop] z;

	//Vector of standard deviations for the shared distributions of the S-R parameters
	vector<lower=0>[2] tau;

	//Hyperparameters of the gaussian process (more detail provided in the model block)
	real<lower=0> rho_dist_MS;
  	real<lower=0> alpha_dist_MS;
  	real<lower=0> sigma_dist_MS;

  	//Standard deviation and correlation of the harvest process error covariance matrix
    real <lower=-0.999,upper=0.999> rho_h;
  	real<lower=0> sigma_h;

  	//Matrix of annual marine survival deviations (from the global mean-reverting random walk process)
  	matrix [n_year-1, n_pop] MS_dev_pop; 
  	matrix [n_year-1, n_pop] harvest_dev; 

  	//Hyperparameters for the hierarchical distribution of hatchery offset terms
  	real mu_offset;
  	vector [n_hatchery] z_offset;
  	real <lower=0> sigma_offset;

  	//Hyperparameters for the hierarchical distribution of initial marine survival terms
  	real mu_init_surv;
	real <lower=0> sigma_init_surv;
	vector[n_pop] z_init_surv;

  	//Hyperparameters for the hierarchical distribution of initial harvest rate terms
	real mu_u_init;
	real <lower=0> sigma_u_init;
	vector[n_pop] z_u_init;

  	//Hyperparameters for the hierarchical distribution of average marine survival terms
	real mu_mu_surv;
	real <lower=0> sigma_mu_surv;
	vector[n_pop] z_mu_surv;

	//Marine survival autocorrelation terms
	vector <lower=-0.99999, upper = 0.99999> [n_pop] phi;

	// Probability of tag recovery for CWT marine survival data (intermediate term for beta-binomial)
	vector <lower=0, upper = 1.0> [n_MS] theta_surv;

    // Smolt and harvest abundance observation error terms
	real <lower=0> sigma_smolt;
	//real <lower=0> sigma_esc;
	real <lower=0> sigma_catch;

}transformed parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  	//Hatchery offset parameter: This parameter is an offset for the likelihood for populations whose marine survival estimates are based off of a nearby hatchery
  	vector[n_hatchery] hatch_offset;

  	//Initial harvest rates for each population
  	vector [n_pop] u_mat_init;   

  	//Annual harvest rates for each population (logit space)
  	matrix [n_year, n_pop] logit_u_mat;    

  	//Annual harvest rates for each population 
	matrix <lower=0, upper = 1> [n_year, n_pop] u_mat;    

  	//Matrix of latent state estimates for smolt abundance by year and population
  	matrix <lower=0> [n_year,n_pop] smolt_est;

	//Matrix of annual early marine survival estimates by population
	matrix [n_year, n_pop] logit_smolt_survival_pop; 

  	//Matrix of latent state estimates for the escapement by year and population
  	matrix <lower=0> [n_year,n_pop] adult_est;

  	//Matrix of latent state estimates for the harvest by year and population
  	matrix <lower=0> [n_year,n_pop] harvest_est;

  	//Matrix of marine survival estimates by year and population
  	matrix <lower=0, upper = 1.0> [n_year,n_pop] smolt_survival;

  	//Matrix of marine survival estimates by year and population
  	matrix <lower=0, upper = 1.0> [n_year,n_pop_MS] smolt_survival_adj;

  	//Standard deviation of recruitment errors for each population
  	vector <lower=0> [n_pop] sigma_R;

  	//Standard deviation of recruitment errors for each population (log transform)
  	vector [n_pop] log_sigma_R;

  	//Log transformation of smolt, escapement, and harvest estimates
  	matrix  [n_year,n_pop] log_smolt_est;
  	matrix  [n_year,n_pop] log_adult_est;
  	matrix  [n_year,n_pop] log_harvest_est;

  	//Productivity term (alpha) and asymptotic recruitment (R_xax) of the Beverton-Holt S-R function by population
  	vector <lower=0> [n_pop] alpha;
  	vector <lower=0> [n_pop] R_max;

  	//Matrix of initial adult abundances for the first two years in which the recruitments arise from unobserved spawning events that are modelled as latent states as in Fleishman et al (2013)
  	matrix [2, n_pop] adult_init;

  	//Correlation matrix for the S-R parameters
  	corr_matrix[2] Omega;

  	//Matrix of the S_R parameters for each population
	matrix[2, n_pop] SR;

	//First year estimate of global (puget MS-wide) marine survival term used to initiate the mean-reverting random walk
  	vector[n_pop] init_surv;	

  	//First year estimate of global (puget MS-wide) marine survival term used to initiate the mean-reverting random walk
  	vector[n_pop] mu_surv;

	//Multivariate non-centered parameterization for drawing the S-R parameters (alpha, beta) from a shared multivariate distribution with mean of gamma, standard deviation tau, 
	//Cholesky factor of the correlation matrix L_omega, and standard scaling factor z
	SR = gamma*u + (diag_pre_multiply(tau,L_Omega) * z);

	//Full correlation matrix calculated by multiplying the cholesky factor of the correlation matrix by its transpose
	Omega = L_Omega * L_Omega';//'

	//Exponentiate from log to normal space the (1) initial adult spawning abundance, (2) alpha S-R parameter, (3) R_max S-R parameter
	adult_init = exp(log_adult_init);
  	alpha = exp(to_vector(SR[1,]));
  	R_max = exp(to_vector(SR[2,])).*stream_dist; // Multiply capacity per KM by stream distance to calculate total capacity

  	//Non-centered hierarchical distribution of recruitment variation term arising from shared distribution among populations
  	for (j in 1:n_pop) log_sigma_R[j] = mu_sigma_R + sigma_sigma_R*z_sigma_R[j];

  	//Non-centered hierarchical distribution of hatchery offset term arising from shared distribution among populations
  	for (j in 1:n_hatchery) hatch_offset[j] = mu_offset + sigma_offset*z_offset[j];

  	//Non-centered hierarchical distribution of initial marine survival terms term arising from shared distribution among populations
  	for (j in 1:n_pop) init_surv[j] = mu_init_surv + sigma_init_surv*z_init_surv[j];
  	
  	//Non-centered hierarchical distribution of initial harvest rates arising from shared distribution among populations
  	for (j in 1:n_pop) u_mat_init[j] = mu_u_init + sigma_u_init*z_u_init[j];

  	//Non-centered hierarchical distribution of average marine survival terms term arising from shared distribution among populations
  	for (j in 1:n_pop) mu_surv[j] = mu_mu_surv + sigma_mu_surv*z_mu_surv[j];

  	//Exponentiate recruitment error standard deviation from log space
  	sigma_R = exp(log_sigma_R);

	//Population-specific marine survival deviations from the global time-series specified as a Gaussian process
	for(i in 1:n_pop){
		for(y in 1:n_year){
			if(y==1){
				logit_smolt_survival_pop[y,i] = init_surv[i];
				logit_u_mat[y,i] = u_mat_init[i];
			}else{
				logit_smolt_survival_pop[y,i] = mu_surv[i] + phi[i]*(logit_smolt_survival_pop[y-1,i]-mu_surv[i]) + MS_dev_pop[y-1, i];
				logit_u_mat[y,i] = logit_u_mat[y-1,i] + harvest_dev[y-1,i];
			}
		}
	}
	
	//Inverse logit transform logit harvest rates
	u_mat = inv_logit(logit_u_mat);

	//Add adjustment term for populations whose marine survival is based on hatchery data
	for (y in 1:n_year){
		for(j in 1:n_hatchery){
			smolt_survival_adj[y,hatchery[j]] = inv_logit(logit_smolt_survival_pop[y,pop_MS[hatchery[j]]] + hatch_offset[j]);
		}
	}

	//Inverse logit transformation of logit smolt survival
	smolt_survival = inv_logit(logit_smolt_survival_pop);

	//no adjustment needed for populations with actual marine survival estimates (not based on hatcheries)
	smolt_survival_adj[,wild]= smolt_survival[,pop_MS[wild]];

	//Beverton-Holt smolt, harvest, and escapement dynamics
	for(i in 1:n_pop){
		for(y in 1:(n_year)){
			if (y< 3){ // First two years of the time series where recruitments arise from unobserved spawning events
				harvest_est[y,i] = adult_init[y, i]*u_mat[y,i]; // harvest calculated as total return times exploitation rate

				adult_est[y,i] = adult_init[y, i]- harvest_est[y,i]; // Escapement calculated as total return minus the harvest

				smolt_est[y,i] = (adult_est[y, i]/(1/alpha[i]  + adult_est[y, i]/R_max[i]))*exp(r_dev[y,i]); // Calculate the resulting smolt recruitment from the escapement according to the Beverton-Holt function

			}else{ //Later years in the time-series where all recruitments arise from observed spawning events

				smolt_est[y,i] = (adult_est[y-2, i]/(1/alpha[i]  + adult_est[y-2, i]/R_max[i]))*exp(r_dev[y,i]); //Smolt recruitment calculated as Beverton-Holt function of spawning abundance 

				harvest_est[y,i] = (smolt_est[y-1,i]*smolt_survival[y-1,i])*u_mat[y,i]; // Harvest is the smolt recruitment, times the marine survival rate of smolts, times the annual explopitation rates

				adult_est[y,i] = smolt_est[y-1,i]*smolt_survival[y-1,i] - harvest_est[y,i]; //Escapement calculated as the smolt recruitment times marine survival minue the harvest
			}
		}
	}

	//Log transform the smolt, escapement, and harvest estimates
	log_smolt_est = log(smolt_est);
	log_harvest_est = log(harvest_est);
	log_adult_est = log(adult_est);

}model {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Transient declaration of (1) Cholesky factor of marine survival covariance matrix, (2) full marine survival covariance matrix, 
//(3) Harvest rate covariance matrix, Choleksy factor of harvest rate covariance matriz
  matrix[n_pop, n_pop] L_K_MS;
  matrix[n_pop, n_pop] K_MS;
  matrix[n_pop, n_pop] K_h;
  matrix[n_pop, n_pop] L_K_h;

// Specify squared exponential kernal for marine survival spatial Gaussian field
  for(i in 1:n_pop){
    for(j in 1:n_pop){
        K_MS[i,j] = square(alpha_dist_MS)*exp(-1*(square(dist[i,j])/(2*square(rho_dist_MS))));
    }
  }

//Add additional error term to the diagonal of the marine survival spatial covariance matrix
  for(n in 1:n_pop){
    K_MS[n,n] = K_MS[n,n] + square(sigma_dist_MS);
  }

  //Cholesky decomposition of marine survival spatial covariance matrix
  L_K_MS = cholesky_decompose(K_MS);

  //Specify covariance matrix for exploitation rate process errors
  for (i in 1:n_pop){
    for(j in 1:n_pop){
        if (i == j){
        K_h[i,j] = sigma_h*sigma_h; //Variance on diagonals
        } else { 
        K_h[i,j] = rho_h*(sigma_h*sigma_h); // Covariance on off-diagonals
      }
    }
  }

  //Cholesky decomposition of harvest rate process errors
  L_K_h = cholesky_decompose(K_h);	

 // Specify multivariate normal distribution of (1) marine survival deviations, (2) harvest rate process errors
  for(y in 1:(n_year-1)){
    MS_dev_pop[y,] ~ multi_normal_cholesky(mu_MS, L_K_MS);
    harvest_dev[y,] ~ multi_normal_cholesky(mu_MS, L_K_h);
  }

   	//Vagua inverse gamma prior on the length scale of the kernel
    rho_dist_MS ~ gamma(1, 0.1);

    // Standard deviation of harvest rate process errors
    sigma_h ~ normal(0,5);

    //Vague normal prior on the marginal SD of the kernel
    alpha_dist_MS ~ normal(0,5);

    //Vague normal prior on the noise term of the GP
    sigma_dist_MS ~ normal(0,5);

   	//Standard normal distribution for the standard scaling factor used in the non-centered hierarchical distrubtion of recruitment deviations
   	z_sigma_R ~ normal(0,1);

   	//Vague normal prior for the mean of the distribution of recruitment standard deviations (sigma_R)
    mu_sigma_R ~ normal(0,5);

    //Vague cauchy prior for the standard deviation of the distribution of recruitment standard deviations (sigma_R)
    sigma_sigma_R ~ normal(0,5);

	for (i in 1:n_pop){
		r_dev[,i] ~ normal(0, sigma_R[i]); // Recruitment deviations lognormally distributed without bias correction factor
		log_adult_init[,i] ~ normal(0, 10); // Vague lognormal prior for initial, unobserved spawning events
	}

	//Scaling factor (z) for the hierarchically/multivariate normally distributed S-R parameters drawn from standard normal prior distribution
	to_vector(z) ~ normal(0,1);

	//Cholesky factor of the correlation matrix drawn from an LKJ prior with shape parameter of 2 (weakly informative in favor of weaker correlation)
	L_Omega ~ lkj_corr_cholesky(2);

	//Priors for the means of the among-population distribution of S-R parameters. Values based on Barrowman et al (2003)
	//alpha
	gamma[1,1] ~ normal(4.27, 2);
	//gamma[1,1] ~ normal(4.27, 5);

	//R_max
	gamma[2,1] ~ normal(7.27, 2);
	//gamma[2,1] ~ normal(7.27, 5);


	//Priors for the standard deviations of the among-population distribution of S-R parameters. Values based on Barrowman et al (2003)
	//alpha
	tau[1] ~ normal(0.43, 0.25);
	//tau[1] ~ normal(0.43, 5);

	//R_max
	tau[2] ~ normal(0.64, 0.25);
	//tau[2] ~ normal(0.64, 5);


	//Initial marine survival hierarchical distribution hyperpriors
	mu_init_surv ~ normal(0,10);
	sigma_init_surv ~ normal(0,5);
	z_init_surv ~ normal(0,1);

	//Average marine survival hierarchical distribution hyperpriors
	mu_u_init ~ normal(0,10);
	sigma_u_init ~ normal(0,5);
	z_u_init ~ normal(0,1);

	//Vague normal prior for the initial year of the global marine survival time series
	mu_mu_surv ~ normal(0,10);
	sigma_mu_surv ~ normal(0,5);
	z_mu_surv ~ normal(0,1);

	//hyperpriors for hatchery offset terms
	z_offset ~ normal(0,1);
	mu_offset ~ normal(0,5);
	sigma_offset ~ normal(0,5);

	//Priors for observation error terms
	sigma_smolt ~ normal(0,5);
	sigma_catch ~ normal(0,5);
	//sigma_esc ~ normal(0,2.5);

	//Likelihoods
	// For each data type, the raw data (appropriately sliced among individual populations from the input vector) is assumed to be distributed 
	//(lognormally in the case of the smolt, harvest, and escapement data, beta-binomially for the CWT marine survival estimates) around the 'true' state estimates, which are produced for each
	// year regardless of whether or not there is data, but will only be evaluated in a likelihood for the years in which such data exists
	for(i in 1:n_pop_smolt){
		smolt_dat[slice_smolt_start[i]:slice_smolt_end[i]] ~ lognormal(log_smolt_est[smolt_true[slice_smolt_start[i]:slice_smolt_end[i]], pop_smolt[i]], sigma_smolt);
	}
	for(i in 1:n_pop_esc){
		esc_dat[slice_esc_start[i]:slice_esc_end[i]] ~ lognormal(log_adult_est[esc_true[slice_esc_start[i]:slice_esc_end[i]], pop_esc[i]], sigma_esc);
	}
	for(i in 1:n_pop_catch){
		harvest_dat[slice_harvest_start[i]:slice_harvest_end[i]] ~ lognormal(log_harvest_est[harvest_true[slice_harvest_start[i]:slice_harvest_end[i]], pop_catch[i]], sigma_catch);
	}
	for(i in 1:n_pop_MS){
		theta_surv[slice_MS_start[i]:slice_MS_end[i]] ~ beta((smolt_survival_adj[MS_true[slice_MS_start[i]:slice_MS_end[i]], i]*(k-2)+1), ((1-smolt_survival_adj[MS_true[slice_MS_start[i]:slice_MS_end[i]],  i])*(k-2)+1));//
		MS_dat_x[slice_MS_start[i]:slice_MS_end[i]] ~ binomial(MS_dat_N[slice_MS_start[i]:slice_MS_end[i]], theta_surv[slice_MS_start[i]:slice_MS_end[i]]);//
	}

} generated quantities {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Declare the hyperparameters (mean and standard deviation) for the the S-R parameters (alpha, and R_max)
	real mu_alpha;
	real sigma_alpha;

	real mu_R_max; 
	real sigma_R_max;

	//Declare individual S-R parameters (alpha, R_max)
	vector [n_pop] log_alpha; 
	vector [n_pop] log_R_max;

	//Declare posterior predictive check parameters for (1) smolt abundance, (2) escapement abundance, (3) harvest abundance
  	matrix  [n_year,n_pop_smolt]  PPC_smolt_abun;
  	matrix  [n_year,n_pop_esc]  PPC_esc_abun;
  	matrix  [n_year,n_pop_catch]  PPC_catch_abun;   
  	vector [n_MS]  PPC_pre_MS;
  	vector [n_MS]  PPC_MS;

  	//Declare model (observation) residuals
  	vector [n_smolt] smolt_resid;   
  	vector [n_esc] esc_resid;
  	vector [n_harvest] harvest_resid;                                    
    vector [n_MS] MS_resid;

    //Vector of adult return forecasts
    vector[n_pop] adult_pred;

    //Spatial correlation matrix
    matrix[n_pop, n_pop] corr_mat;

    //Calculate correlation by distance (as opposed to covariance), useful for easier plotting
    for(i in 1:n_pop){
    	for(j in 1:n_pop){
       		corr_mat[i,j] = exp(-1*(square(dist[i,j])/(2*square(rho_dist_MS))));
    	}
  	}

    //Forecast of next year's adult return
    for(i in 1:n_pop){
		adult_pred[i] = smolt_est[n_year,i]*smolt_survival[n_year,i];
	}

	//Pull the standard deviations for the hyperdistributions of alpha and R_max from the tau vector
	sigma_alpha = tau[1];
	sigma_R_max = tau[2];

	//Pull the means of the hyperdistributions of alpha and R_max from the gamma matrix
	mu_alpha = gamma[1,1];
	mu_R_max = gamma[2,1];

	// Pull the individual S-R parameters (alpha, R_max) from the S-R matrix used for the Multivariate normal/Cholesky implementation of the Beverton-Holt Parameters
	log_alpha = to_vector(SR[1,]);
	log_R_max = to_vector(SR[2,]);

	//Compute posterior predictive distributions for (1) smolt abundance, (2) escapement abundance, (3) harvest abundance, and (4) marine survival
	for(i in 1:n_pop_smolt){
		for(y in 1:n_year){
			PPC_smolt_abun[y,i] = lognormal_rng(log_smolt_est[y,pop_smolt[i]], sigma_smolt);
		}
	}
	for(i in 1:n_pop_esc){
		for(y in 1:n_year){
		    PPC_esc_abun[y,i] = lognormal_rng(log_adult_est[y,pop_esc[i]], sigma_esc);
		 }
	}
	for(i in 1:n_pop_catch){
		for(y in 1:n_year){
			PPC_catch_abun[y,i] = lognormal_rng(log_harvest_est[y,pop_catch[i]], sigma_catch);

		}
	}
	//Posterior predictive check for the beta-binomial CWT data must be specified differently because without the sample size, we can only generate predictive distributions for years in which the data exist
	for(i in 1:n_pop_MS){
		for(j in slice_MS_start[i]:slice_MS_end[i]){
		   PPC_pre_MS[j] = beta_rng((smolt_survival_adj[MS_true[j], i]*(k-2)+1), ((1-smolt_survival_adj[MS_true[j], i])*(k-2)+1));
		   PPC_MS[j] = binomial_rng(MS_dat_N[j], PPC_pre_MS[j]);
		}
	}

	//Calculate model (observation residuals)
	for(i in 1:n_pop_smolt){
		smolt_resid[slice_smolt_start[i]:slice_smolt_end[i]] = to_vector(smolt_dat[slice_smolt_start[i]:slice_smolt_end[i]]) - smolt_est[smolt_true[slice_smolt_start[i]:slice_smolt_end[i]], pop_smolt[i]];
	}
	for(i in 1:n_pop_esc){
		esc_resid[slice_esc_start[i]:slice_esc_end[i]] = to_vector(esc_dat[slice_esc_start[i]:slice_esc_end[i]]) - adult_est[esc_true[slice_esc_start[i]:slice_esc_end[i]], pop_esc[i]];
	}
	for(i in 1:n_pop_catch){
		harvest_resid[slice_harvest_start[i]:slice_harvest_end[i]] = to_vector(harvest_dat[slice_harvest_start[i]:slice_harvest_end[i]]) - harvest_est[harvest_true[slice_harvest_start[i]:slice_harvest_end[i]], pop_catch[i]];
	}
	for(i in 1:n_pop_MS){
		MS_resid[slice_MS_start[i]:slice_MS_end[i]] = to_vector(to_vector(MS_dat_x[slice_MS_start[i]:slice_MS_end[i]])./to_vector(MS_dat_N[slice_MS_start[i]:slice_MS_end[i]])) - (smolt_survival_adj[MS_true[slice_MS_start[i]:slice_MS_end[i]], i]);
	}
}
