[ PROB ]
hmctdmr :: phenytoin, one-compartment with Michaelis-Menten elimination (matches inst/stan/phenytoin.stan)

[ SET ]
req = "Vmax Norm_Vmax Km VD CL IPRED AUC", end = 30, rtol = 1e-5, delta = 0.01

[ GLOBAL ]
#define CP (CENT/VD)

[ CMT ] @annotated
CENT : Central compartment (mg)
AUC  : AUC

[ PARAM ] @annotated
WT       : 60   : Body weight (kg)
LBW      : 60   : Lean body weight (kg)
CrCL     : 1    : Creatinine clearance (mL/min)
Vmax     : 500  : Maximum elimination rate (mg/day) for 70 kg
Km       : 5    : Michaelis constant (mg/L)
VD_NR    : 0.8  : Volume (L/kg LBW)
CL_SLOPE : 0.01 : Renal clearance slope, fixed as in the Stan model

[ MAIN ]
double Norm_Vmax = Vmax * pow(WT / 70, 0.6);
double VD = VD_NR * LBW;
double CL = CL_SLOPE * CrCL * 60 / 1000;

[ ODE ]
dxdt_CENT = -(CL + (Norm_Vmax / 24) / (Km + CP)) * CP;
dxdt_AUC  = CP;
if (SS_ADVANCE) dxdt_AUC = 0;

[ TABLE ]
double IPRED = CP;

[ CAPTURE ] @annotated
Vmax      : Maximum elimination rate (mg/day) for 70 kg
Norm_Vmax : Maximum elimination rate (mg/day) for the patient, Vmax * (WT/70)^0.6
Km    : Michaelis constant (mg/L)
VD    : Volume (L)
CL    : Linear clearance (L/h)
IPRED : Individual-predicted concentration (mg/L)
