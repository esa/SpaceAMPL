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
	InitialPositionx2: x2[1] = x02;
	InitialPositiony2: y2[1] = y02;
	InitialPositionz2: z2[1] = z02;
	FinalPositionx2  : x2[n2] = xf2;
	FinalPositiony2  : y2[n2] = yf2;
	FinalPositionz2  : z2[n2] = zf2;
#-----------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f1{i in I2} = dx2[i];
var f2{i in I2} = dy2[i];
var f3{i in I2} = dz2[i];
var f4{i in I2} = -x2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2);
var f5{i in I2} = -y2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2);
var f6{i in I2} = -z2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm{i in J2} = (x2[i]+x2[i+1])/2 + tf2/(n2-1)/8 * (f1[i] - f1[i+1]);
var ym{i in J2} = (y2[i]+y2[i+1])/2 + tf2/(n2-1)/8 * (f2[i] - f2[i+1]);
var zm{i in J2} = (z2[i]+z2[i+1])/2 + tf2/(n2-1)/8 * (f3[i] - f3[i+1]);
var dxm{i in J2} = (dx2[i]+dx2[i+1])/2 + tf2/(n2-1)/8 * (f4[i] - f4[i+1]);
var dym{i in J2} = (dy2[i]+dy2[i+1])/2 + tf2/(n2-1)/8 * (f5[i] - f5[i+1]);
var dzm{i in J2} = (dz2[i]+dz2[i+1])/2 + tf2/(n2-1)/8 * (f6[i] - f6[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m{i in J2} = dxm[i];
var f2m{i in J2} = dym[i];
var f3m{i in J2} = dzm[i];
var f4m{i in J2} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f5m{i in J2} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f6m{i in J2} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx2{i in J2}: x2[i] = x2[i+1] - tf2/(n2-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicy2{i in J2}: y2[i] = y2[i+1] - tf2/(n2-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicz2{i in J2}: z2[i] = z2[i+1] - tf2/(n2-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicdx2{i in J2}: dx2[i] = dx2[i+1] - tf2/(n2-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamicdy2{i in J2}: dy2[i] = dy2[i+1] - tf2/(n2-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicdz2{i in J2}: dz2[i] = dz2[i+1] - tf2/(n2-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
		
#--------------------------------------------------------------------------
#Initial Guess for the LP
let {i in I2} x2[i] :=    x02*cos(tf2/(n2-1)*(i-1)) - y02*sin(tf2/(n2-1)*(i-1));
let {i in I2} y2[i] :=    x02*sin(tf2/(n2-1)*(i-1)) + y02*cos(tf2/(n2-1)*(i-1));
let {i in I2} dx2[i] :=  -x02*sin(tf2/(n2-1)*(i-1)) - y02*cos(tf2/(n2-1)*(i-1));
let {i in I2} dy2[i] :=   x02*cos(tf2/(n2-1)*(i-1)) - y02*sin(tf2/(n2-1)*(i-1));
#-----------------------------------------------------------------------

solve;

#--------------------------------------------------------------------------
#Initial Guess for the remaining variables
let {i in I2} m2[i] := 1;
#The following guess is needed for snopt as otherwise snopt starts with an
#arithmetic error
let {i in I2} ux2[i]:=Tmax*0.0001;
let {i in I2} uy2[i]:=Tmax*0.0001;
let {i in I2} uz2[i]:=Tmax*0.0001;
let {i in J2} uxm2[i]:=Tmax*0.0001;
let {i in J2} uym2[i]:=Tmax*0.0001;
let {i in J2} uzm2[i]:=Tmax*0.0001;

#-----------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the Initial Guess
for {i in I2} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e,%17.16e\n",x2[i],y2[i],z2[i],dx2[i],dy2[i],dz2[i],m2[i],ux2[i],uy2[i],uz2[i]>out/InitialGuess2.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n",  ti2/d2u , tF2-ti2/d2u > out/TimesGuess2.out;

#--------------------------------------------------------------------------
#As to calculate the LP initial guess we made use of a simplified dynamic we 
#we need to dispose some variables (note that AMPL will automatically dispose
#also those variables linked to the one listed here


purge
f1,f2,f3,f4,f5,f6,f1m,f2m,f3m,f4m,f5m,f6m,InitialPositionx2,InitialPositiony2,InitialPositionz2,FinalPositionx2,FinalPositiony2,FinalPositionz2;
#--------------------------------------------------------------------------

