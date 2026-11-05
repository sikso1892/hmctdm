functions {
  vector OneCptIV_Model(real t, vector cmt, array[] real parms, array[] real x_r, array[] int x_i){
    
    real CL = parms[1];
    real VD = parms[2];
    real ke     = CL/VD;

    vector[1] dxdt_cmt;
    
    dxdt_cmt[1] = -ke * cmt[1];

    return dxdt_cmt;
  }

  real cv_to_sd(real cv) {
    return sqrt(log(cv^2+1));
  }

  real sigma(real conc){
    return (conc * 0.15 + 0.25);
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
  
  array[nt] real<lower=0> CrCL;
  array[nt] real<lower=0> LBW;

  real Prior_CL_NR;
  real Prior_VD_NR; 
  real Prior_CL_SLOPE;
 
  real Prior_CL_NR_omega;
  real Prior_VD_NR_omega;
  real Prior_CL_SLOPE_omega;
  
}

transformed data{
  int nTheta = 2;
  int nCmt = 1;
}

parameters {

  real<lower=0> CL_NR;
  real<lower=0> VD_NR;
  real<lower=0> CL_SLOPE;

}

transformed parameters {
  array[nTheta] real theta;

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;
  array[nObs] real<lower=0> sigmaEPS;
  
  real CL = (CL_NR *  max(LBW)  + CL_SLOPE * max(CrCL) ) *  60 / 1000;
  real VD = (VD_NR *  max(LBW));
  
  theta[1] = CL;
  theta[2] = VD;
    
  matrix[nCmt, nt] x = pmx_solve_rk45(OneCptIV_Model, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[1, ] ./ VD;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

  CL_NR ~ lognormal(log(Prior_CL_NR), Prior_CL_NR_omega);
  VD_NR ~ lognormal(log(Prior_VD_NR), Prior_VD_NR_omega);
  CL_SLOPE ~ lognormal(log(Prior_CL_SLOPE), Prior_CL_SLOPE_omega);

  for ( i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
  }
}
