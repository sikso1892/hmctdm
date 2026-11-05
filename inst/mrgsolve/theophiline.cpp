[ PROB ]
hmctdmr :: theophiline, one-compartment oral (matches inst/stan/theophiline.stan)

[ SET ]
req = "CL VC IPRED AUC", end = 30, rtol = 1e-5, delta = 0.5

[ GLOBAL ]
#define CP (CENT/VC)

[ CMT ] @annotated
GUT  : Absorption compartment (mg)
CENT : Central compartment (mg)
AUC  : AUC

[ PARAM ] @annotated
WT    : 70   : Body weight (kg)
LBW   : 70   : Lean body weight (kg)
CrCL  : 30   : Creatinine clearance (mL/min)
ka    : 0.27 : Absorption rate constant (1/h), fixed as in the Stan model
CL_NR : 40   : Clearance (mL/h/kg LBW)
VC_NR : 0.5  : Volume (L/kg LBW)

[ MAIN ]
double CL = CL_NR * LBW / 1000;
double VC = VC_NR * LBW;
double ke = CL / VC;

[ ODE ]
dxdt_GUT  = -ka * GUT;
dxdt_CENT = ka * GUT - ke * CENT;
dxdt_AUC  = CP;
if (SS_ADVANCE) dxdt_AUC = 0;

[ TABLE ]
double IPRED = CP;

[ CAPTURE ] @annotated
CL    : Clearance (L/h)
VC    : Volume (L)
IPRED : Individual-predicted concentration (mg/L)
