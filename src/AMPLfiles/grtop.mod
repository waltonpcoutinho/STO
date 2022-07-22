#----------------------------------------------------------------------------------------------#
#---------------------------------------"ADD DESCRIPTION"--------------------------------------#
#----------------------------------------------------------------------------------------------#
reset;

#----------------------------------------------------------------------------------------------#
#----------------------------------------"GET FILE NAME"---------------------------------------#
#----------------------------------------------------------------------------------------------#
#"pause command"
#shell 'read -rs';

#"sets up an output file name"
#param outputFile = ("Solutions/" & instanceName) symbolic;

#----------------------------------------------------------------------------------------------#
#---------------------------------"FINISH READ DATA FROM FILE"---------------------------------#
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
#-------------------------------------"DEFINE PARAMETERS"--------------------------------------#
#----------------------------------------------------------------------------------------------#
param N;
param delta = 1e-3;

#"define time related constants"
set Time := {0..N-1} ordered;
set States := {1..6};
set Controls := {1..2};

#"define the system's matrices"
param A{States,States};
param B{States,Controls};
param normA;
param normB;

#"equilibrium conditions"
param Yeq{States};
param Ueq{Controls};

#"initial and final conditions"
param Y0{States};
param U0{Controls};
param Yf{States};
param Uf{Controls};

#"lower and upper bounds"
param Ylb{States};
param Ulb{Controls};
param Yub{States};
param Uub{Controls};

#"upper bound on norm of derivatives"
param norm_yDotub;
param norm_uDotub;

#"auxiliary variable for computing the norm on eps"
param overall_norm = normA*norm_yDotub + normB*norm_uDotub;

#"lower bound on flight time"
param tflb;

#"vertices locations"
param xBar;
param yBar;
param hBar;
param rBar;
param coneMinZ;
param coneMaxZ;

#"landing zones locations"
param xTilde;
param yTilde;
param hTilde;
param rTilde;

#"photographing parameters"
param UbarH;
param barH;
param hatGamma;
param hatMu;

#----------------------------------------------------------------------------------------------#
#--------------------------------------"DEFINE VARIABLES"--------------------------------------#
#----------------------------------------------------------------------------------------------#
#"state variables"
var Y{Time,States};
var U{Time,Controls};

#"discretization error"
var epsilon{Time,States};
var z{Time} >= 1e-3;

#"final time and step size"
var tf >= 0;
var h = tf/(N-1);

#----------------------------------------------------------------------------------------------#
#-------------------------------------"WRITE OPT. PROBLEM"-------------------------------------#
#----------------------------------------------------------------------------------------------#

minimize objective: tf + sum{t in Time} z[t]; 

subject to

#"dynamic Euler constraints"
EOMs{t in 0..(N-2), s in States}:
Y[t+1,s] = Y[t,s] + h*(sum{j in States} A[s,j]*(Y[t,j] - Yeq[j])
                   + sum{k in Controls}  B[s,k]*(U[t,k] - Ueq[k])) + epsilon[t,s];

#"initial states"
InitY{s in States}:  Y[0,s] = Y0[s];

#"initial controls"
InitC{c in Controls}:  U[0,c] = U0[c];

#"bounds on state variables"
boundsY{t in Time, s in States}:     Ylb[s] <= Y[t,s] <= Yub[s];

#"bounds on control variables"
boundsC{t in Time, c in Controls}:   Ulb[c] <= U[t,c] <= Uub[c];

#"constraints on the norm of epsilon"
normEps1{t in Time, s in States}: z[t] >=   epsilon[t,s];
normEps2{t in Time, s in States}: z[t] >= - epsilon[t,s];

#"upper bound on the local truncation errors"
auxEpsBounds{t in Time}: z[t] <= 0.5*(h^2)*overall_norm;

#"lower bound on the flight time"
tf_lower_bound: tf >= tflb;

#visiting and landing constraints
#in separate files

#"stabilisation constraints"
perturbation1{t in Time, c in Controls}: 
    U[t,c] >= (t/(N-1))*(Ueq[c] - delta) + Ulb[c]*(1 - (t/(N-1)));
perturbation2{t in Time, c in Controls}: 
    U[t,c] <= (t/(N-1))*(Ueq[c] + delta) + Uub[c]*(1 - (t/(N-1)));

#perturbation{t in 1..N-1, c in Controls}: 
#   Ueq[c] - delta <= U[t,c] <= Ueq[c] + delta; 






