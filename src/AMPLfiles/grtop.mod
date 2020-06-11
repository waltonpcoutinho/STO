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
param xi = 0.001;

#"define time related constants"
set Time := {0..N-1} ordered;
set States := {1..6};
set Controls := {1..2};

#"define the system's matrices"
param A{States,States};
param B{States,Controls};

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

#"vertices locations"
param xBar;
param yBar;
param hBar;
param rBar;
param coneMinZ;
param coneMaxZ;

#landing zones locations
param xTilde;
param yTilde;
param hTilde;
param rTilde;

#photographing parameters
param UbarH;
param barH;
param hatGamma;
param hatMu;

#----------------------------------------------------------------------------------------------#
#--------------------------------------"DEFINE VARIABLES"--------------------------------------#
#----------------------------------------------------------------------------------------------#
#state variables
var Y{Time,States};
var U{Time,Controls};

#discretization error
var epsilon >= 0;

#final time and step size
var tf >= 0;
var h = tf/(N-1);

#----------------------------------------------------------------------------------------------#
#-------------------------------------"WRITE OPT. PROBLEM"-------------------------------------#
#----------------------------------------------------------------------------------------------#
minimize objective: epsilon; 

subject to

#"dynamic Euler constraints"
EOMs1{t in 0..(N-2), s in States}:
Y[t+1,s] >= Y[t,s] + h*(sum{j in States} A[s,j]*(Y[t,j] - Yeq[j])
                   + sum{k in Controls}  B[s,k]*(U[t,k] - Ueq[k])) - epsilon;
EOMs2{t in 0..(N-2), s in States}:
Y[t+1,s] <= Y[t,s] + h*(sum{j in States} A[s,j]*(Y[t,j] - Yeq[j])
                   + sum{k in Controls}  B[s,k]*(U[t,k] - Ueq[k])) + epsilon;

#"initial states"
InitY{s in States}:  Y[0,s] = Y0[s];

#"initial controls"
InitC{c in Controls}:  U[0,c] = U0[c];

#"bounds on state variables"
boundsY{t in Time, s in States}:     Ylb[s] <= Y[t,s] <= Yub[s];

#"bounds on control variables"
boundsC{t in Time, c in Controls}:   Ulb[c] <= U[t,c] <= Uub[c];

#"stabilisation constraints"
perturbation1{t in Time, c in Controls}: 
    U[t,c] >= (t/(N-1))*(Ueq[c] - xi) + Ulb[c]*(((N-1)-t)/(N-1));
perturbation2{t in Time, c in Controls}: 
    U[t,c] <= (t/(N-1))*(Ueq[c] + xi) + Uub[c]*(((N-1)-t)/(N-1));








