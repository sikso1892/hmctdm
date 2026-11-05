functions {

  vector ODE_NonLinear(real t, vector cmt, array[] real parms, array[] real x_r, array[] int x_i){
    
    real Vmax = parms[1];
    real Km = parms[2];
    real VD = parms[3];
    real CL = parms[4];

    vector[1] dxdt_cmt;

    real CP = cmt[1]/VD;
    
    dxdt_cmt[1] = -( CL + (Vmax/24) / (Km + CP) ) * CP;

    return dxdt_cmt;
  }
  real cv_to_sd(real cv) {
    return sqrt(log(cv^2+1));
  }

  real sigma(real conc){
    return (conc * 0.1 + 1.0);
  }

}

data {

  int<lower=1> nt;                        // number of events
  int<lower=1> nObs;
  array[nObs] int<lower=1> iObs;
  vector<lower=0>[nObs] cObs;
  
  // torstan parameters (nonmem data)
  array[nt] int<lower=1> cmt;

  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;

  array[nt] real CrCL;
  array[nt] real LBW;
  array[nt] real WT;

  real Prior_Vmax;
  real Prior_Km;
  real Prior_VD_NR; 
 
  real Prior_Vmax_omega;
  real Prior_Km_omega;
  real Prior_VD_NR_omega;
}

transformed data{

  int nCmt = 1;
  real CL_SLOPE=0.01;

}
parameters{
  real ln_Vmax;
  real ln_Km;
  real ln_VD_NR;

}
transformed parameters  {

  array[4] real parms;

  real Vmax = exp(ln_Vmax);
  real Km = exp(ln_Km);
  real VD_NR = exp(ln_VD_NR);
  
  real CL = CL_SLOPE * CrCL[1] * 60 / 1000;
  real VD = VD_NR * LBW[1];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;

  array[nObs] real<lower=0> sigmaEPS;
  
  parms[1] = Vmax * pow(WT[1]/70, 0.6);
  parms[2] = Km;
  parms[3] = VD;
  parms[4] = CL;

  matrix[nCmt, nt] x = pmx_solve_rk45(ODE_NonLinear, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, parms, 1e-5, 1e-8, 1e5);   
  
  cHat = x[1, ] ./ VD;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

   /* ... declarations ... statements ... */
   ln_Vmax ~ normal(log(Prior_Vmax), Prior_Vmax_omega);
   ln_Km ~ normal(log(Prior_Km), Prior_Km_omega);
   ln_VD_NR ~ normal(log(Prior_VD_NR), Prior_VD_NR_omega);

   for (i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
   }

}
