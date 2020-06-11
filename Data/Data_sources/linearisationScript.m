#!/usr/local/bin/MathematicaScript -script
(*how to run this script:
	math -noprompt -run '<<linearisationScript.m'*)

(*Include packages*)
<< ToMatlab.m;
Needs["MatrixManipulation`"];

(*ODE *)
dX = f1 = v0*Cos[gamma0]*Sin[phi0] + betaC*z0;
dY = f2 = v0*Cos[gamma0]*Cos[phi0];
dZ = f3 = v0*Sin[gamma0];  
dXdd = f4 = -Sin[gamma0]*Sin[phi0]*v0*dGamma + Cos[gamma0]*Cos[phi0]*v0*dPhi + Cos[gamma0]*Sin[phi0]*dV;
dYdd = f5 = -Cos[phi0]*Sin[gamma0]*v0*dGamma - Cos[gamma0]*Sin[phi0]*v0*dPhi + Cos[gamma0]*Cos[phi0]*dV;
dZdd = f6 = Cos[gamma0]*v0*dGamma + Sin[gamma0]*dV;
dV = f7 = -((rho*S (CD0 + K*CL0^2)*v0^2)/(2 mass)) - gravity*Sin[gamma0] - betaC*v0*Sin[gamma0]*Cos[gamma0]*Sin[phi0];
dGamma = f8 = ((rho*S*CL0*v0*Cos[mu0])/(2*mass)) - gravity*Cos[gamma0]/v0 + betaC*Sin[phi0]*Sin[gamma0]^2;
dPhi = f9 = - ((rho*S*CL0*v0^2*Sin[mu0])/(2*mass*v0*Cos[gamma0])) - (betaC*Sin[gamma0] Cos[phi0]/Cos[gamma0]); 

(*system of ODEs*)
M = {{f1}, {f2}, {f3}, {f4}, {f5}, {f6}, {f7}, {f8}, {f9}};
(*compute ODEs 1st derivatives*)
col1 = D[M, x0];
col2 = D[M, y0];
col3 = D[M, z0];
col4 = D[M, dX0];
col5 = D[M, dY0];
col6 = D[M, dZ0];
col7 = D[M, v0];
col8 = D[M, gamma0];
col9 = D[M, phi0];
col10 = D[M, CL0];
col11 = D[M, mu0]; 
jacobian = 
  AppendRows[col1, col2, col3, col4, col5, col6, col7, col8, col9, 
   col10, col11] ;
Dimensions[jacobian];
(*create state and control matrices and export them*)

matrixA = jacobian[[1 ;; 9, 1 ;; 9]];
Dimensions[matrixA];
matrixB = jacobian[[1 ;; 9, 10 ;; 11]];
Dimensions[matrixB];
exp1 = ToMatlab[matrixA, "exp1", 1000000]; 
exp2 = ToMatlab[matrixB, "exp2", 1000000];
file = OpenWrite[
   "/home/walton/Dropbox/Univesity_of_Southampton/Doutorado/Problems/\
GliderRouting/Implementation/Mathematica/ODE_linear/tomatlab_test1.m"];
WriteMatlab[exp1, file];
WriteMatlab[exp2, file];   
Close[file];

(*solve system of ODEs to get equilibrium points*)

data = ReadList[
   "/home/walton/Dropbox/Univesity_of_Southampton/Doutorado/Problems/\
GliderRouting/Implementation/LinearisedModel/Data/GliderData.dat"];
rho = data[[1]];
CD0 = data[[2]];
K = data[[3]];
gravity = data[[4]];
mass = data[[5]];
S = data[[6]];
betaC = data[[9]];
(*stationary points*)  	
(*practical*)

eq = ReadList[
   "/home/walton/Dropbox/Univesity_of_Southampton/Doutorado/Problems/\
GliderRouting/Implementation/LinearisedModel/Data/EqPoints.dat"];
v0 = eq[[1]];
gamma0 = eq[[2]];
phi0 = eq[[3]];
mu0 = eq[[4]];
CL0 := Sqrt[CD0/K];

(*output State and Control Matrices*)

Export["/home/walton/Dropbox/Univesity_of_Southampton/Doutorado/\
Problems/GliderRouting/Implementation/LinearisedModel/Data/\
StateMatrix.dat", matrixA, "TSV"];
Export["/home/walton/Dropbox/Univesity_of_Southampton/Doutorado/\
Problems/GliderRouting/Implementation/LinearisedModel/Data/\
ControlMatrix.dat", matrixB, "TSV"];
matrixA // MatrixForm;
matrixB // MatrixForm;

(*Clear all variables*)
ClearAll["Global`*"];
Quit[];




