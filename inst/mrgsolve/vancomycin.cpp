[ PROB ]
hmctdmr :: vancomycin, two-compartment IV (matches inst/stan/vancomycin.stan)

[ SET ]
req = "CL VC VP IPRED AUC", end = 30, rtol = 1e-5, delta = 0.5

[ GLOBAL ]
#define CP (CENT/VC)

[ CMT ] @annotated
CENT   : Central compartment (mg)
PERIPH : Peripheral compartment (mg)
AUC    : AUC

[ PARAM ] @annotated
WT       : 70      : Body weight (kg)
LBW      : 70      : Lean body weight (kg)
CrCL     : 30      : Creatinine clearance (mL/min)
CL_NR    : 0.05    : Non-renal clearance (mL/min/kg LBW)
CL_SLOPE : 0.75    : Renal clearance slope
VC_NR    : 0.21    : Central volume (L/kg LBW)
k12      : 1.12    : Central -> peripheral rate constant (1/h)
k21      : 0.48    : Peripheral -> central rate constant (1/h)
UNIT_CL  : 0.06    : mL/min -> L/h

[ MAIN ]
double CL = (CL_NR * LBW + CL_SLOPE * CrCL) * UNIT_CL;
double VC = VC_NR * LBW;
double VP = k12 / k21 * VC;
double ke = CL / VC;

[ ODE ]
dxdt_CENT   = k21 * PERIPH - (ke + k12) * CENT;
dxdt_PERIPH = k12 * CENT - k21 * PERIPH;
dxdt_AUC    = CP;
if (SS_ADVANCE) dxdt_AUC = 0;

[ TABLE ]
double IPRED = CP;

[ CAPTURE ] @annotated
CL    : Clearance (L/h)
VC    : Central volume (L)
VP    : Peripheral volume (L)
IPRED : Individual-predicted concentration (mg/L)
