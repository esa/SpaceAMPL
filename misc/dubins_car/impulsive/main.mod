#Optimal steering for the Dubins car
#In this problem transcription the analytical solution of the ballistic problem is 
#used as to obtain a reduced dimensional problem as in the MAILTO algorithm
#Written by Dario Izzo (February 2009)
#

#Parameters---------------------------------------------------
param n:=100;			#Numbers of nodes
#-------------------------------------------------------------

#Sets---------------------------------------------------------
set I := 1..n;
set J := 1..n-1;
#-------------------------------------------------------------

#Variables----------------------------------------------------
var tf, >=0;					#final time
var x{i in I};			#state 1
var y{i in I};			#state 2
var theta{i in I};		#state 3
var Dtheta{i in J};		#control 1
#-------------------------------------------------------------

var dt = tf/(n-1);		#

#Objective----------------------------------------------------
	minimize time: tf;
#-------------------------------------------------------------

#Constraints--------------------------------------------------

	#Dynamic explicit solution constraints
	subject to 
		dynamic1{i in J}: x[i+1] = x[i] + dt * cos(theta[i] + Dtheta[i]);
		dynamic2{i in J}: y[i+1] = y[i] + dt * sin(theta[i] + Dtheta[i]);
		dynamic3{i in J}: theta[i+1] = theta[i] + Dtheta[i];
		
	#Boundary Conditions
	subject to InitialPositionx: x[1] = 10;
	subject to InitialPositiony: y[1] = 0;
	subject to InitialPositiontheta: theta[1] = 2;
	subject to FinalPositionx: x[n] = 0;
	subject to FinalPositiony: y[n] = 0;
	subject to FinalPositiontheta: theta[n]=2;

	#Conrol Constraint
	subject to controlmagnitude1{i in J}:  Dtheta[i] - dt <=0;
	subject to controlmagnitude2{i in J}:  Dtheta[i] + dt >=0;
#-------------------------------------------------------------

#-------------------------------------------------------------	
#Initial Conditions
	let tf := 12;
	let {i in I} theta[i]:=3.14;
#-------------------------------------------------------------

#-------------------------------------------------------------
#Solver Options
option solver snopt;
option substout 1;
option show_stats 1;
options snopt_options "outlev=2 Major_iterations=1500";
#-------------------------------------------------------------

#-------------------------------------------------------------
#Solve!!!
solve;
#display {i in I} y[i,1];
#display {i in I} y[i,2];
#display {i in I} y[i,3];
#display {i in I} u[i];
#display {i in J} um[i];
#-------------------------------------------------------------
	
	
	
	
	
	
