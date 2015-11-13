#Single phase transfer (iniguess.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to build an initial guess for the transfer



#--------------------------------------------------------------------------
#Problem definition
subject to 
	InitialPositionx1: x1[1] = x01;
	InitialPositiony1: y1[1] = y01;
	InitialPositionz1: z1[1] = z01;
	FinalPositionx1  : x1[n1] = xf1;
	FinalPositiony1  : y1[n1] = yf1;
	FinalPositionz1  : z1[n1] = zf1;
#-----------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f1{i in I1} = dx1[i];
var f2{i in I1} = dy1[i];
var f3{i in I1} = dz1[i];
var f4{i in I1} = -x1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2);
var f5{i in I1} = -y1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2);
var f6{i in I1} = -z1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm{i in J1} = (x1[i]+x1[i+1])/2 + tf1/(n1-1)/8 * (f1[i] - f1[i+1]);
var ym{i in J1} = (y1[i]+y1[i+1])/2 + tf1/(n1-1)/8 * (f2[i] - f2[i+1]);
var zm{i in J1} = (z1[i]+z1[i+1])/2 + tf1/(n1-1)/8 * (f3[i] - f3[i+1]);
var dxm{i in J1} = (dx1[i]+dx1[i+1])/2 + tf1/(n1-1)/8 * (f4[i] - f4[i+1]);
var dym{i in J1} = (dy1[i]+dy1[i+1])/2 + tf1/(n1-1)/8 * (f5[i] - f5[i+1]);
var dzm{i in J1} = (dz1[i]+dz1[i+1])/2 + tf1/(n1-1)/8 * (f6[i] - f6[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m{i in J1} = dxm[i];
var f2m{i in J1} = dym[i];
var f3m{i in J1} = dzm[i];
var f4m{i in J1} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f5m{i in J1} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f6m{i in J1} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx1{i in J1}: x1[i] = x1[i+1] - tf1/(n1-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicy1{i in J1}: y1[i] = y1[i+1] - tf1/(n1-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicz1{i in J1}: z1[i] = z1[i+1] - tf1/(n1-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicdx1{i in J1}: dx1[i] = dx1[i+1] - tf1/(n1-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamicdy1{i in J1}: dy1[i] = dy1[i+1] - tf1/(n1-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicdz1{i in J1}: dz1[i] = dz1[i+1] - tf1/(n1-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
		
#--------------------------------------------------------------------------
#Initial Guess for the LP
let {i in I1} x1[i] :=    x01*cos(tf1/(n1-1)*(i-1)) - y01*sin(tf1/(n1-1)*(i-1));
let {i in I1} y1[i] :=    x01*sin(tf1/(n1-1)*(i-1)) + y01*cos(tf1/(n1-1)*(i-1));
let {i in I1} dx1[i] :=  -x01*sin(tf1/(n1-1)*(i-1)) - y01*cos(tf1/(n1-1)*(i-1));
let {i in I1} dy1[i] :=   x01*cos(tf1/(n1-1)*(i-1)) - y01*sin(tf1/(n1-1)*(i-1));
#-----------------------------------------------------------------------

solve;

#--------------------------------------------------------------------------
#Initial Guess for the remaining variables
let {i in I1} m1[i] := 1;

#The following guess is needed for snopt as otherwise snopt starts with an
#arithmetic error
let {i in I1} ux1[i]:=Tmax*0.0001;
let {i in I1} uy1[i]:=Tmax*0.0001;
let {i in I1} uz1[i]:=Tmax*0.0001;
let {i in J1} uxm1[i]:=Tmax*0.0001;
let {i in J1} uym1[i]:=Tmax*0.0001;
let {i in J1} uzm1[i]:=Tmax*0.0001;

#-----------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the Initial Guess
for {i in I1} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e,%17.16e\n",x1[i],y1[i],z1[i],dx1[i],dy1[i],dz1[i],m1[i],ux1[i],uy1[i],uz1[i]>out/InitialGuess1.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti1/d2u , tF1-ti1/d2u > out/TimesGuess1.out;



#--------------------------------------------------------------------------
#As to calculate the LP initial guess we made use of a simplified dynamic we 
#we need to dispose some variables (note that AMPL will automatically dispose
#also those variables linked to the one listed here


purge
f1,f2,f3,f4,f5,f6,f1m,f2m,f3m,f4m,f5m,f6m,InitialPositionx1,InitialPositiony1,InitialPositionz1,FinalPositionx1,FinalPositiony1,FinalPositionz1;
#--------------------------------------------------------------------------

