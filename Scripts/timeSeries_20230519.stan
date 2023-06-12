
data {
  int<lower=0> Ncreek;
  int<lower=0> Ntime;
  int<lower=0> Nobs; //total observations
  int<lower=0> Nstations; //probably 2, downstream and upstream
  int<lower=0> Nspecies;
  int<lower=0> Nconstruction; //N culvert states; removed (=2) vs. nonremoved (=1)
  int<lower=0> MinconstructionTimepoint; //first timepoint with construction
  int<lower=0> NconstructionTimepoints;
  int<lower=0> time_idx[Nobs];
  int<lower=0> creek_idx[Nobs];
  int<lower=0> species_idx[Nobs];
  int<lower=0> station_idx[Nobs]; //downstream = 1, upstream = 2
  int<lower=0> construction_idx[Nobs];
  vector[Nobs] y_logeDNA;
  
  // int<lower=0>N_unobserved;
  // int<lower=0>unobserved_time_idx[N_unobserved];
  // int<lower=0>unobserved_creek_idx[N_unobserved];
  // int<lower=0>unobserved_station_idx[N_unobserved];
  // int<lower=0>unobserved_species_idx[N_unobserved];

}

parameters {
  
  //real mu_0[Nstations, Nspecies, Ncreek]; //starting concentration for eDNA (log scale)
  vector[Ntime] alpha[Ncreek, Nspecies]; //creek effect
  real epsilon[Ntime, Nstations, Nspecies, Ncreek]; //autocorrelation effect
  vector<lower=-1, upper = 1>[Nspecies] beta_1; //slope relating present concentration to past concentration
  real<lower=0> phi; //process variance term on autoregression
  
  real eta_unk[Ntime, Nspecies, Ncreek]; //culvert effect; true parameter portion
  real mu_eta[Nspecies]; //hier; mean of culvert effects
  // real<lower=0> sigma_eta;  //hier; variance among culvert effects
  // real<lower=0> sigma_gamma; //hier; variance among construction effects
  
  vector[Nspecies] mu_alpha; //hier; mean of creek effects
  real<lower=0> sigma_alpha; //hier; variance among creek effects
  
  vector<lower=0>[Nspecies] sigma_dna; //observation variance
}


transformed parameters {
  real mu[Ntime, Nstations, Nspecies, Ncreek]; //mean eDNA (log scale) from which obs was drawn // dimensions = c(time, station, creek) 
  real eta[Ntime, Nstations, Nspecies, Ncreek]; // culvert effect
  
  /////// eta (culvert effect) as mix of parameters and fixed values, where downstream site/times are defined as 0.0
      eta[,1,,] = rep_array(0.0, Ntime, Nspecies, Ncreek);
      //then fill in the rest of the array w param values
       for (t in 1:Ntime){
          for (j in 1:Nspecies){
            for (i in 1:Ncreek){
              eta[t,2,j,i] = eta_unk[t,j,i];
            }}}
            
  ///////Now, calculate expected values (mu)
    for (t in 1:Ntime){
      for (d in 1:Nstations){
        for (j in 1:Nspecies){
          for (i in 1:Ncreek){
            mu[t,d,j,i] =  
                alpha[i, j, t] + // creek effect
                eta[t, d, j, i] + // upstream-downstream effect; defined as 0 where d == 1
                epsilon[t, d, j, i];  //autocorrelation; depends upon mu at previous timestep
          }
        }
      }
    }
}


model {

  for (i in 1:Nobs){
    y_logeDNA[i] ~ normal(mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
                          sigma_dna[species_idx[i]]); 
  }


    //Sample autoregression effect, with mean  
    for (t in 2:Ntime){
      for (d in 1:Nstations){
        for (j in 1:Nspecies){
            for (i in 1:Ncreek){                                  
              if(i>1){
                epsilon[t,d,j,i] ~ normal(beta_1[j] * mu[t-1,d,j,i], phi);  //  note: object mu has one more timestep than object epsilon; the first index of mu reflects the first timestep; the first index of epsilon reflects the second timestep.  
              }
              if(i==1){
                epsilon[t,d,j,i] ~ normal(0, pow(pow(phi,2) / (1-pow(beta_1[j],2)),-2));  // 
                // This is the asymptotic variance for the the first time step.
                // note: object mu has one more timestep than object epsilon; the first index of mu reflects the first timestep; the first index of epsilon reflects the second timestep.
              }
            }
          }
        }}



  // Priors
  // for (d in 1:Nstations){
  //   for (i in 1:Ncreek){
  //     for (j in 1:Nspecies){
  //       epsilon[,d,j,i] ~ normal(0, 5); //   // //keep positive?
  //       mu_0[d,j,i] ~ normal(0, 5);
  //     }}}
  
      beta_1 ~ uniform(-1,1);  //autocorr slope
      
      for (i in 1:Ncreek){
      for (j in 1:Nspecies){
        alpha[i,j] ~ normal(mu_alpha[j], sigma_alpha); //creek effect
        eta_unk[,j,i] ~ normal(mu_eta[j], 1); //culvert effect
        // gamma_unk[,j] ~ normal(0, 2); //construction effect
      }
      }
      
    mu_eta ~ normal(0,5);
    // sigma_eta ~ gamma(1,2); //var culvert effect
    // sigma_gamma ~ gamma(1,2); //var construction effect
    mu_alpha ~ normal(0,10); //mean of creek effects; hier
    sigma_alpha ~ gamma(1,2); //var creek effect
    phi ~ gamma(1,2); //var autocorr effect
    sigma_dna ~ gamma(1,1); //var observations
  
}

generated quantities{
//  real delta[Ntime-1, Nspecies, Ncreek]; //difference between upstream and downstream of culvert
// 
// // for (i in 1:Nobs){
// //   delta[i] = normal_rng(mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
// //                           sigma_dna[species_idx[i]]);
// // }
// 
//     for (t in 2:Ntime){
//         for (j in 1:Nspecies){
//           for (i in 1:Ncreek){
//              for (r in 1:Nconstruction){
//               delta[t-1,j,i] = eta[t-1,2,j,i] - eta[t-1,1,j,i]; //upstream minus downstream
//             }
//           }
//         }
//       }



 // vector[Nobs] log_lik;
 //   for (i in 1:Nobs){
 //    log_lik[i]  =  normal_lpdf(y_logeDNA[i] | mu[time_idx[i], station_idx[i], species_idx[i], creek_idx[i]], 
 //                          sigma_dna[species_idx[i]]);
 // 
 //  }



 
}
