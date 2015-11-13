#--------------------------------------------------------------------------
#Single phase transfer (ocp.inc)
#Written by Dario Izzo (February 2009)
#
#This file may be included in main.mod and generates and solves the full 
#optimal control problem
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Problem definition. Only the boundary constraint and thrust magnitude 
#define entirely the problem
subject to 
	ControlMagnitude{i in 2..n-1} : uT[i]  <= Tmax * tf/(n-1);
	FinalC3:	 VINFf <= sqrt(C3arr);
	InitialC3:	 VINF <= sqrt(C3dep);
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
minimize 
	obj1: sum{i in 2..n-1} uT[i]^2;	#quadratic control
solve;
drop obj1;

maximize
	obj2: exp(-(sum{i in 2..n-1} uT[i]));	#mass
solve;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Print the solution x,y,z,dx,dy,dz,m,ux,uy,uz variables and more

var m{i in 2..n-1} = exp(-(sum{j in 2..i} uT[j])/Isp/g0);
include include/writesolution.inc;

#Print the initial and final times
printf "%17.16e  %17.16e \n", ti/d2u , tF-ti/d2u > out/Times.out;
close out/Times.out;
#Print the final mass
printf "%17.16e \n", m[n-1] > out/mass.out;
close out/mass.out;
#--------------------------------------------------------------------------

