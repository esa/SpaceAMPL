#Optimal steering for the Dubins car
#Trapezoidal Formula is used to set the dynamic constraints
#
#Written by Dario Izzo (October 2007)
#

#Parameters---------------------------------------------------
param n:=50;			#Numbers of nodes
#-------------------------------------------------------------

#Sets---------------------------------------------------------
set I := 1..n;
set J := 1..n-1;
#-------------------------------------------------------------

#Variables----------------------------------------------------
var tf;
var x{i in I};
var y{i in I};
var theta{i in I};
var u{i in I};
#-------------------------------------------------------------

#Constraints--------------------------------------------------

	#Boundary Conditions
	subject to InitialPositionx: x[1] = 10;
	subject to InitialPositiony: y[1] = 0;
	subject to InitialPositiontheta: theta[1] = 2;
	subject to FinalPositionx: x[n] = 0;
	subject to FinalPositiony: y[n] = 0;
	subject to FinalPositionTheta: theta[n] = 2;

	#Control Constraint
	subject to controlmagnitude{i in 1..n}: -1 <= u[i] <= 1;
	
	var dt = tf / (n-1);

	#Dynamical Constraints (trapezoidal collocation)
	subject to dynamicx{i in J}: 
		 x[i+1] = x[i] + dt * (cos(theta[i]) + cos(theta[i+1])) / 2;
	subject to dynamicy{i in J}: 
		 y[i+1] = y[i] + dt * (sin(theta[i]) + sin(theta[i+1])) / 2;
	subject to dynamictheta{i in J}:
		 theta[i+1] = theta[i] + dt * (u[i] + u[i+1]) / 2;

#-------------------------------------------------------------	

#Initial Conditions-------------------------------------------
	let tf := 12;
	let {i in I} theta[i]:=3.14;
#-------------------------------------------------------------

#Objective----------------------------------------------------
	minimize time: tf;
#-------------------------------------------------------------

#-------------------------------------------------------------
#Solver Options
option solver snopt;
option substout 1;
option show_stats 1;
options snopt_options "outlev=2 Major_iterations=1500";
#-------------------------------------------------------------

#Solve!!!
solve;
#display _varname, _var;
