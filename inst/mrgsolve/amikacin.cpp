[ PROB ]
hmctdmr :: amikacin, one-compartment IV (matches inst/stan/amikacin.stan)

[ SET ]
req = "CL VD IPRED AUC", end = 30, rtol = 1e-5, delta = 0.5

[ GLOBAL ]
#define CP (CENT/VD)

[ CMT ] @annotated
CENT : Central compartment (mg)
AUC  : AUC

[ PARAM ] @annotated
WT       : 70      : Body weight (kg)
LBW      : 70      : Lean body weight (kg)
CrCL     : 30      : Creatinine clearance (mL/min)
CL_NR    : 0.0417  : Non-renal clearance (mL/min/kg LBW)
CL_SLOPE : 0.815   : Renal clearance slope
VD_NR    : 0.27    : Volume (L/kg LBW)
UNIT_CL  : 0.06    : mL/min -> L/h

[ MAIN ]
double CL = (CL_NR * LBW + CL_SLOPE * CrCL) * UNIT_CL;
double VD = VD_NR * LBW;

[ ODE ]
dxdt_CENT = -CL * CP;
dxdt_AUC  = CP;
if (SS_ADVANCE) dxdt_AUC = 0;

[ TABLE ]
double IPRED = CP;

[ CAPTURE ] @annotated
CL    : Clearance (L/h)
VD    : Volume (L)
IPRED : Individual-predicted concentration (mg/L)
