#Optimal steering for the Dubins car
#Hermite-Simpson Formula is used to set the dynamic constraints
#Matrix notation is used to write the dynamic, this produces 
#a larger problem as presolve and substitute are not able to
#get rid of the extra variables introduced
#Written by Dario Izzo (October 2007)
#

#Parameters---------------------------------------------------
param n:=50;			#Numbers of nodes
param N:=3;			#Number of equations
#-------------------------------------------------------------

#Sets---------------------------------------------------------
set I := 1..n;
set J := 1..n-1;
set Neq := 1..N;
#-------------------------------------------------------------

#Variables----------------------------------------------------
var tf;				#final time
var y{i in I,k in Neq};		#state
var u{i in I};			#controls
var ym{i in J, k in Neq};	#state in midpoints
var um{i in J};			#controls in midpoint	
var f{i in I, k in Neq};	#dynamic
var fm{i in J, k in Neq};	#dynamic in midpoints
#-------------------------------------------------------------

#Objective----------------------------------------------------
	minimize time: tf;
#-------------------------------------------------------------

#Constraints--------------------------------------------------

	#Dynamic definition constraints (presolve should substitute most of these)
	subject to 
		dynamic1{i in I}: f[i,1] = cos(y[i,3]);
		dynamic2{i in I}: f[i,2] = sin(y[i,3]);
		dynamic3{i in I}: f[i,3] = u[i];
	
	subject to 
		HermiteInterpolation{i in J, j in Neq}: 
			ym[i,j] = (y[i,j]+y[i+1,j])/2 + tf/(n-1)/8 * (f[i,j]-f[i+1,j]);
	
	subject to 
		dynamicm1{i in J}: fm[i,1] = cos(ym[i,3]);
		dynamicm2{i in J}: fm[i,2] = sin(ym[i,3]);
		dynamicm3{i in J}: fm[i,3] = um[i];


	#Simpson formula constraints
	subject to SimpsonRule{i in J, j in Neq}:
		y[i,j] = y[i+1,j] - tf/(n-1)/6 * (f[i,j] + f[i+1,j] + 4*fm[i,j]);
		
	#Boundary Conditions
	subject to InitialPositionx: y[1,1] = 10;
	subject to InitialPositiony: y[1,2] = 0;
	subject to InitialPositiontheta: y[1,3] = 2;
	subject to FinalPositionx: y[n,1] = 0;
	subject to FinalPositiony: y[n,2] = 0;
	subject to FinalPositiontheta: y[n,3]=2;

	#Conrol Constraint
	subject to controlmagnitude{i in I}: -1 <= u[i] <= 1;
	subject to controlmagnitudem{i in J}: -1 <= um[i] <= 1;

#-------------------------------------------------------------

#-------------------------------------------------------------	
#Initial Conditions
	let tf := 12;
	let {i in I} y[i,3]:=3.14;
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
	
	
	
	
	
	
