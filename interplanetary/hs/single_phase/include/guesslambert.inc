#Single phase transfer (iniguess.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to build an initial guess for the transfer

fix timod;
fix tfmod;

#--------------------------------------------------------------------------
#Problem definition
subject to 
	InitialPositionx: x[1] = x0;
	InitialPositiony: y[1] = y0;
	InitialPositionz: z[1] = z0;
	FinalPositionx  : x[n] = xf;
	FinalPositiony  : y[n] = yf;
	FinalPositionz  : z[n] = zf;
#-----------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f1{i in I} = dx[i];
var f2{i in I} = dy[i];
var f3{i in I} = dz[i];
var f4{i in I} = -x[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2);
var f5{i in I} = -y[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2);
var f6{i in I} = -z[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm{i in J} = (x[i]+x[i+1])/2 + tf/(n-1)/8 * (f1[i] - f1[i+1]);
var ym{i in J} = (y[i]+y[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
var zm{i in J} = (z[i]+z[i+1])/2 + tf/(n-1)/8 * (f3[i] - f3[i+1]);
var dxm{i in J} = (dx[i]+dx[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
var dym{i in J} = (dy[i]+dy[i+1])/2 + tf/(n-1)/8 * (f5[i] - f5[i+1]);
var dzm{i in J} = (dz[i]+dz[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m{i in J} = dxm[i];
var f2m{i in J} = dym[i];
var f3m{i in J} = dzm[i];
var f4m{i in J} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f5m{i in J} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
var f6m{i in J} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx{i in J}: x[i] = x[i+1] - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicy{i in J}: y[i] = y[i+1] - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicz{i in J}: z[i] = z[i+1] - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicdx{i in J}: dx[i] = dx[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamicdy{i in J}: dy[i] = dy[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicdz{i in J}: dz[i] = dz[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
		
#--------------------------------------------------------------------------
#Initial Guess for the LP
let {i in I} x[i] :=   x0*cos(tf/(n-1)*(i-1)) - y0*sin(tf/(n-1)*(i-1));
let {i in I} y[i] :=   x0*sin(tf/(n-1)*(i-1)) + y0*cos(tf/(n-1)*(i-1));
let {i in I} dx[i] :=  -x0*sin(tf/(n-1)*(i-1)) - y0*cos(tf/(n-1)*(i-1));
let {i in I} dy[i] :=  x0*cos(tf/(n-1)*(i-1)) - y0*sin(tf/(n-1)*(i-1));
#-----------------------------------------------------------------------

solve;

#--------------------------------------------------------------------------
#Initial Guess for the remaining variables
let {i in I} m[i] := 1;

#--------------------------------------------------------------------------
#The following guess is needed for snopt as otherwise snopt starts with an
#arithmetic error
let {i in I} ux[i]:=Tmax*0.00001;
let {i in I} uy[i]:=Tmax*0.00001;
let {i in I} uz[i]:=Tmax*0.00001;
let {i in J} uxm[i]:=Tmax*0.00001;
let {i in J} uym[i]:=Tmax*0.00001;
let {i in J} uzm[i]:=Tmax*0.00001;
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the Initial Guess
for {i in I} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",x[i],y[i],z[i],dx[i],dy[i],dz[i],m[i],ux[i],uy[i],uz[i] >out/InitialGuess.out;
}
close out/InitialGuess.out;
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Print the initial and final times
printf "%17.16e, %17.16e \n", ti/d2u , tF-ti/d2u > out/TimesGuess.out;
close out/TimesGuess.out;
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#As to calculate the LP initial guess we made use of a simplified dynamic we 
#we need to dispose some variables (note that AMPL will automatically dispose
#also those variables linked to the one listed here

unfix timod;
unfix tfmod;
purge f1,f2,f3,f4,f5,f6,f1m,f2m,f3m,f4m,f5m,f6m,InitialPositionx,InitialPositiony,InitialPositionz,FinalPositionx,FinalPositiony,FinalPositionz;
#--------------------------------------------------------------------------

