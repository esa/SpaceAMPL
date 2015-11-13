#A simple Lambert problem transcribed as a constraint satisfaction problem
#The initial guess determines the multiple revolution selected
#Written by Dario Izzo (November 2006)

#Parameters-------------------------------------------------------------
param n:=20;                			#Number of points
#-----------------------------------------------------------------------

#Sets-------------------------------------------------------------------
set I := 1..n;
set J := 1..n-1;
#-----------------------------------------------------------------------

#Variables--------------------------------------------------------------
var x1{i in I};
var x2{i in I};
var x3{i in I};
var x4{i in I};
var x5{i in I};
var x6{i in I};
#-----------------------------------------------------------------------


#Lambert problem definition---------------------------------------------
	param tf:=20;
	subject to 
		InitialPositionx: x1[1] = 1;
		InitialPositiony: x2[1] = 0;
		InitialPositionz: x3[1] = 0;
		FinalPositionx  : x1[n] = 0;
		FinalPositiony  : x2[n] = 1;
		FinalPositionz  : x3[n] = 0;
#-----------------------------------------------------------------------		

#Dynamic definition at grid points--------------------------------------
var f1{i in I} = x4[i];
var f2{i in I} = x5[i];
var f3{i in I} = x6[i];
var f4{i in I} = -x1[i]/(x1[i]**2+x2[i]**2+x3[i]**2)^(3/2);
var f5{i in I} = -x2[i]/(x1[i]**2+x2[i]**2+x3[i]**2)^(3/2);
var f6{i in I} = -x3[i]/(x1[i]**2+x2[i]**2+x3[i]**2)^(3/2);
#-----------------------------------------------------------------------

#State definition at mid-points via Simpson interpolation---------------
var x1m{i in J} = (x1[i]+x1[i+1])/2 + tf/(n-1)/8 * (f1[i] - f1[i+1]);
var x2m{i in J} = (x2[i]+x2[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
var x3m{i in J} = (x3[i]+x3[i+1])/2 + tf/(n-1)/8 * (f3[i] - f3[i+1]);
var x4m{i in J} = (x4[i]+x4[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
var x5m{i in J} = (x5[i]+x5[i+1])/2 + tf/(n-1)/8 * (f5[i] - f5[i+1]);
var x6m{i in J} = (x6[i]+x6[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
#-----------------------------------------------------------------------

#Dynamic definition at mid-points---------------------------------------
var f1m{i in J} = x4m[i];
var f2m{i in J} = x5m[i];
var f3m{i in J} = x6m[i];
var f4m{i in J} = -x1m[i]/(x1m[i]**2+x2m[i]**2+x3m[i]**2)^(3/2);
var f5m{i in J} = -x2m[i]/(x1m[i]**2+x2m[i]**2+x3m[i]**2)^(3/2);
var f6m{i in J} = -x3m[i]/(x1m[i]**2+x2m[i]**2+x3m[i]**2)^(3/2);
#-----------------------------------------------------------------------

#Hermite Formula--------------------------------------------------------

	subject to 
		dynamicx1{i in J}: x1[i] = x1[i+1] - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
		dynamicx2{i in J}: x2[i] = x2[i+1] - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
		dynamicx3{i in J}: x3[i] = x3[i+1] - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
		dynamicx4{i in J}: x4[i] = x4[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
		dynamicx5{i in J}: x5[i] = x5[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
		dynamicx6{i in J}: x6[i] = x6[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
		

#Initial Guess----------------------------------------------------------
let {i in I} x1[i] := cos(tf/(n-1)*(i-1));
let {i in I} x2[i] := sin(tf/(n-1)*(i-1));
let {i in I} x4[i] := -sin(tf/(n-1)*(i-1));
let {i in I} x5[i] := cos(tf/(n-1)*(i-1));
#-----------------------------------------------------------------------

#Solver Options---------------------------------------------------------
option solver snopt;
#option substout 1;		#Substitutes redundant variables reducing the dimension
option show_stats 1;		#Prints the result of the substitution and of presolve
#option linelim 1;
options snopt_options "outlev=2 Major_iterations=1500";
#------------------------------------------------------------------------

#Solve!!!---------------------------------------------------------------
solve;
#-----------------------------------------------------------------------

#Print the solution-----------------------------------------------------
#The State
for {i in I} {
printf "%e, %e, %e, %e, %e, %e \n",x1[i],x2[i],x3[i],x4[i],x5[i],x6[i] > State.out;
}

#------------------------------------------------------------------------
